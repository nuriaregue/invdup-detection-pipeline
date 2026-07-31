"""
1. Per cada read guardar
    - Si overlap: segments on es produeix overlap
    - Si GAP: reg1 -SAs - reg2

2. Buscar events 
Per cada read de INV_DUP mirar
    - Si overlap: mirar si final de la regió 1 i inici de la regió 2 són similars 
    TEMA DELS STRANDS: clar realment inici i fi estan invertits en les regions -
        vale segons l'strand es mira quina és la coordenada d'inici o fi
    - Si GAP: mirar si reg1 final, SAs i regió 2 inici son similars
Similars seria unes coordenades amb +=Marge i = chr
La cosa es com es sabria amb quines coordenades mirar si son similars? 
"""
from statistics import median
import pandas as pd
from pathlib import Path
from scripts.config import (
    MIN_GEOMETRY_MAPQ,
    MIN_SA_RECOVERY_MAPQ,
    COORDINATE_MARGIN
)

from scripts.reads_analysis.classify_geometry import find_sa_supported_region_gap


def process_candidate_coordinates(candidates_tsv:Path,segments_tsv:Path,pairs_tsv:Path):
    if candidates_tsv.exists() and segments_tsv.exists() and pairs_tsv.exists():
        bam_name = candidates_tsv.parent.name

        candidates_df=pd.read_csv(candidates_tsv,sep="\t")
        segments_df=pd.read_csv(segments_tsv,sep="\t")
        pairs_df=pd.read_csv(pairs_tsv,sep="\t")
        candidates_inv_dup=candidates_df[candidates_df["CANDIDATE_TYPE"] == "POSSIBLE_INV_DUP"]
        rows = []
        signatures={}

        for candidate in candidates_inv_dup.itertuples(index=False):
            bam = candidate.BAM
            inversion = candidate.INVERSION
            read_id = candidate.READ_ID
            classification_reason = candidate.CLASSIFICATION_REASON
            read_pairs = pairs_df[(pairs_df["BAM"] == bam)&(pairs_df["READ_ID"] == read_id)]

            signatures[read_id]={
                "bam": bam,
                "inversion": inversion,
                "classification_reason": classification_reason,
                "type": "UNKNOWN"
            }
            read_segments = segments_df[(segments_df["BAM"] == bam)&(segments_df["READ_ID"] == read_id)]
            if classification_reason in ["COMPATIBLE_INTRACHROMOSOMAL_OVERLAP","COMPATIBLE_TWO_SEGMENT_PATTERN"]:
                best_pair = get_overlap_support(read_pairs)
                if best_pair is not None:
                    signatures[read_id]["type"]="LOCAL_OVERLAP"
                    segment_a,segment_b=get_pair_segments(best_pair,read_segments)
                    create_local_overlap_entry(segment_a,segment_b,signatures[read_id])
                    for segment in [segment_a,segment_b]: #segment_a and segment_b are dicts

                        rows.append({
                            "BAM": bam,
                            "INVERSION": inversion,
                            "READ_ID": read_id,
                            "CLASSIFICATION_REASON": classification_reason,
                            "ROLE": "OVERLAP_SUPPORT",
                            "ALIGNMENT_INDEX": segment["ALIGNMENT_INDEX"],
                            "SOURCE": segment["SOURCE"],
                            "STRAND": segment["STRAND"],
                            "REF_CHR": segment["REF_CHR"],
                            "REF_START": segment["REF_START"],
                            "REF_END": segment["REF_END"]
                        })

            #Once we have the REGION segments we add SUPLEMENTARY alignments that explain the gaps
            elif classification_reason in ["DISTANT_INTRACHROMOSOMAL_GAP_WITH_INVERTED_SA","INTERCHROMOSOMAL_GAP_WITH_INVERTED_SA","REGION_GAP_WITH_INVERTED_SA_SUPPORT"]:
                region_1,region_2,sa_support=get_sa_gap_support(read_segments)
                if region_1 is not None:
                    signatures[read_id]["type"]="GAP"
                    create_gap_flanks_entry(region_1,region_2,signatures[read_id])
                    for region in [region_1,region_2]: #region_1 and region_2 are dicts
                        rows.append({
                            "BAM": bam,
                            "INVERSION": inversion,
                            "READ_ID": read_id,
                            "CLASSIFICATION_REASON": classification_reason,
                            "ROLE": "REGION_FLANK",
                            "ALIGNMENT_INDEX": region["ALIGNMENT_INDEX"],
                            "SOURCE": region["SOURCE"],
                            "STRAND": region["STRAND"],
                            "REF_CHR": region["REF_CHR"],
                            "REF_START": region["REF_START"],
                            "REF_END": region["REF_END"]
                        })
                    signatures[read_id]["sa_support"]=[]
                    for sa in sa_support:
                        create_gap_sa_entry(sa,signatures[read_id]["sa_support"])
                        rows.append({
                            "BAM": bam,
                            "INVERSION": inversion,
                            "READ_ID": read_id,
                            "CLASSIFICATION_REASON": classification_reason,
                            "ROLE": "SA_GAP_SUPPORT",
                            "ALIGNMENT_INDEX": sa["ALIGNMENT_INDEX"],
                            "SOURCE": sa["SOURCE"],
                            "STRAND": sa["STRAND"],
                            "REF_CHR": sa["REF_CHR"],
                            "REF_START": sa["REF_START"],
                            "REF_END": sa["REF_END"]
                        })

        overlap_events = create_overlap_events(signatures)
        gap_events =create_gap_events(signatures)
        events_df = summarize_coordinate_events(signatures,overlap_events,gap_events)
        events_path = candidates_tsv.parent / f"{bam_name}_coordinate_events.tsv"
        events_df.to_csv(events_path,sep="\t",index=False)
        return signatures, events_df
    else:
        raise FileNotFoundError(f"Missing input file in {candidates_tsv.parent}")

def create_local_overlap_entry(segment_a:dict,segment_b:dict,read_dict:dict):
    """
    Method to create a local_overlap entry to the signatures dictionary
    """
    ref_start_a=segment_a["REF_START"]
    ref_start_b=segment_b["REF_START"]
    ref_end_a=segment_a["REF_END"]
    ref_end_b=segment_b["REF_END"]
    overlap_start = max(ref_start_a, ref_start_b)
    overlap_end = min(ref_end_a, ref_end_b)
    chrm=segment_a["REF_CHR"]
    read_dict["chr"]=chrm
    read_dict["overlap_start"]=overlap_start
    read_dict["overlap_end"]=overlap_end
    
def create_gap_flanks_entry(region_1,region_2,read_dict:dict):
    """
    Method to create a gap entry (flanks) to the signatures dictionary
    """
    strand = region_1["STRAND"]
    chrm=region_1["REF_CHR"]
    if strand=="-":
        flank_1=region_1["REF_START"]
        flank_2=region_2["REF_END"]
    else:
        flank_1=region_1["REF_END"]
        flank_2=region_2["REF_START"]
    
    read_dict["region_start"] = min(flank_1, flank_2)
    read_dict["region_end"] = max(flank_1, flank_2)
    read_dict["chr"]=chrm
    read_dict["strand"]=strand
def create_gap_sa_entry(sa:dict,sa_list:list):
    """
    Method to create a gap entry (sa) to the signatures dictionary
    """
    sa_dict={}
    sa_dict["alignment_index"]=sa["ALIGNMENT_INDEX"]
    sa_dict["sa_start"]=sa["REF_START"]
    sa_dict["sa_end"]=sa["REF_END"]
    sa_dict["chr"]=sa["REF_CHR"]
    sa_dict["strand"]=sa["STRAND"]
    sa_list.append(sa_dict)

def get_overlap_support(read_pairs):
    """
    Method to select the overlapping segment pair with the strongest reference overlap support.
    """
    compatible_pairs = read_pairs[(read_pairs["SAME_CHROMOSOME"] == True)
        & (read_pairs["OPPOSITE_STRANDS"] == True)
        & (read_pairs["REFERENCE_OVERLAP_BP"] > 0)]
    if compatible_pairs.empty:
        return None
    #If there are more than one overlapping pair. the larger one is chosen with idxmax
    best_pair = compatible_pairs.loc[
        compatible_pairs["REFERENCE_OVERLAP_BP"].idxmax()
    ]
    return best_pair

def get_sa_gap_support(read_segments):
    #TODO: arreglar el refactoring
    """
    Method to extract the regions flanks and supplementary segments that support a query gap.
    """
    geometry_alignments = read_segments[read_segments["MAPQ"] >= MIN_GEOMETRY_MAPQ]

    supplementary_alignments = read_segments[read_segments["SOURCE"].isin(["RECOVERED_SA","SA_TAG_ONLY"])
        & (read_segments["MAPQ"] >= MIN_SA_RECOVERY_MAPQ)]
    region_1,region_2,sa_support,gap_type = find_sa_supported_region_gap(geometry_alignments,supplementary_alignments)

    return region_1,region_2,sa_support
    
def get_pair_segments(best_pair: pd.Series,read_segments: pd.DataFrame):
    """
    Method to extract the two read segments corresponding to the selected overlap pair.
    """
    index_a = best_pair["ALIGNMENT_INDEX_A"]
    index_b = best_pair["ALIGNMENT_INDEX_B"]
    segment_a = read_segments[read_segments["ALIGNMENT_INDEX"] == index_a].iloc[0].to_dict() #We convert it to dict to have the same format that sa_regions
    segment_b = read_segments[read_segments["ALIGNMENT_INDEX"] == index_b].iloc[0].to_dict() 

    return segment_a, segment_b

def choose_event_anchor(signatures:list,start_name:str,end_name:str):
    """
    Method to select the most recurrent overlap coordinate.
    To use it with regions and overlaps the start end end key name is introduced as a parameter
    """
    start_matches = 0
    end_matches = 0

    for i in range(len(signatures)):
        for j in range(i + 1, len(signatures)):
            signature_1 = signatures[i] 
            signature_2 = signatures[j]

            if abs(signature_1[start_name] - signature_2[start_name]) <= COORDINATE_MARGIN:
                start_matches += 1

            if abs(signature_1[end_name] - signature_2[end_name]) <= COORDINATE_MARGIN:
                end_matches += 1
    
    if start_matches == 0 and end_matches == 0:
        return None
    elif start_matches >= end_matches:
        return start_name
    else:
        return end_name

def create_overlap_events(signatures:dict):
    """
    Method to group local overlap signatures into recurrent coordinate events.
    """
    overlap_reads = get_signatures_by_type(signatures,"LOCAL_OVERLAP")

    recurring_anchor_type = choose_event_anchor(list(overlap_reads.values()),"overlap_start","overlap_end")

    events = {}

    if recurring_anchor_type is not None:
        sorted_reads = sort_signatures(overlap_reads,recurring_anchor_type)
        current_event_id = None
        
        for read_id, signature in sorted_reads:
            recurring_coord = signature[recurring_anchor_type]
            overlap_size = signature["overlap_end"] - signature["overlap_start"]
            same_event = False

            if current_event_id is not None:
                current_event = events[current_event_id]
                same_chr = signature["chr"] == current_event["chr"]
                close_coordinate = abs(recurring_coord - current_event["anchor"]) <= COORDINATE_MARGIN
                same_event = same_chr and close_coordinate
            if same_event:
                current_event["reads"].append(read_id)
                current_event["anchor_values"].append(recurring_coord)
                current_event["overlap_sizes"].append(overlap_size)
                current_event["anchor"] = median(current_event["anchor_values"])

            else:
                current_event_id = f"OVERLAP_EVENT_{len(events) + 1}"

                events[current_event_id] = {
                    "inversion": signature["inversion"],
                    "chr": signature["chr"],
                    "anchor_type": recurring_anchor_type,
                    "anchor": recurring_coord,
                    "anchor_values": [recurring_coord],
                    "overlap_sizes": [overlap_size],
                    "reads": [read_id]
                }

    return events
def is_compatible_overlap_event(signature:dict,event:dict,recurring_coord):
    """
    Method to check if a LOCAL_OVERLAP signature is compatible with an existing event.
    """
    same_chr = signature["chr"] == event["chr"]
    close_coordinate = abs(recurring_coord - event["anchor"]) <= COORDINATE_MARGIN

    return same_chr and close_coordinate

def have_similar_sa(signature_1:dict, signature_2:dict):
    """
    Method to check if two GAP signatures have similar SA support.
    """
    sa_support_1 = signature_1["sa_support"]
    sa_support_2 = signature_2["sa_support"]
    #GAP signatures need supplementary support in both reads to be comparable.
    if len(sa_support_1) == 0 or len(sa_support_2) == 0:
        return False
    #We require all SA alignments from the read with fewer supplementary segments to find a compatible match.
    required_matches = min(len(sa_support_1),len(sa_support_2))

    #We keep a separate list of available SA alignments so the same SA cannot be matched more than once.
    matched_sa = 0
    available_sa = sa_support_2.copy()
    for sa_1 in sa_support_1:
        matched_sa_2 = None
        for sa_2 in available_sa:
            if matched_sa_2 is None:
                same_chr = sa_1["chr"] == sa_2["chr"]
                close_start = abs(sa_1["sa_start"] - sa_2["sa_start"]) <= COORDINATE_MARGIN
                close_end = abs(sa_1["sa_end"] - sa_2["sa_end"]) <= COORDINATE_MARGIN
                #We consider two SA alignments similar if at least one reference boundary is recurrent.
                if same_chr and (close_start or close_end):
                    matched_sa_2 = sa_2

        if matched_sa_2 is not None:
            matched_sa += 1
            available_sa.remove(matched_sa_2) #The matched SA is removed so it cannot support another SA from the same read.

    return matched_sa == required_matches

def create_gap_events(signatures:dict):
    """
    Method to group GAP signatures into recurrent coordinate events.
    """
    gap_reads = get_signatures_by_type(signatures,"GAP")
    recurring_anchor_type = choose_event_anchor(list(gap_reads.values()), "region_start","region_end")

    events = {}

    if recurring_anchor_type is not None:
        sorted_reads = sort_signatures(gap_reads,recurring_anchor_type)
        for read_id, signature in sorted_reads:
            add_gap_to_event(read_id,signature,events,recurring_anchor_type)

    return events

def get_signatures_by_type(signatures:dict,signature_type:str):
    """
    Method to select signatures of one type.
    """
    selected = {}

    for read_id, signature in signatures.items():
        if signature["type"] == signature_type:
            selected[read_id] = signature

    return selected

def sort_signatures(signatures:dict,anchor_type:str):
    """
    Method to sort signatures by anchor coordinate.
    """
    sorted_signatures = sorted(signatures.items(),
        key=lambda read: (read[1][anchor_type])
    )

    return sorted_signatures

def add_gap_to_event(read_id:str,signature:dict,events:dict,anchor_type:str):
    """
    Method to add a GAP signature to a compatible event or create a new one.
    """
    recurring_coord = signature[anchor_type]
    compatible_event = None

    for event_id, event in events.items():
        if compatible_event is None:
            if is_compatible_gap_event(signature,event,recurring_coord):
                compatible_event = event_id

    if compatible_event is not None:
        event = events[compatible_event]

        event["reads"].append(read_id)
        event["anchor_values"].append(recurring_coord)
        event["anchor"] = median(event["anchor_values"])

    else:
        event_id = f"GAP_EVENT_{len(events) + 1}"

        events[event_id] = {
            "inversion": signature["inversion"],
            "chr": signature["chr"],
            "strand": signature["strand"],
            "anchor_type": anchor_type,
            "anchor": recurring_coord,
            "anchor_values": [recurring_coord],
            "reads": [read_id],
            "representative_signature": signature
        }

def is_compatible_gap_event(signature:dict,event:dict,recurring_coord):
    """
    Method to check if a GAP signature is compatible with an existing event.
    """
    close_coordinate = abs(recurring_coord - event["anchor"]) <= COORDINATE_MARGIN
    similar_sa = have_similar_sa(signature,event["representative_signature"])

    return close_coordinate and similar_sa

def summarize_coordinate_events(signatures:dict, overlap_events:dict, gap_events:dict):
    """
    Method to summarize the coordinate events detected in one BAM.
    """
    rows = []

    for event_id, event in overlap_events.items():
        event_signatures = []

        for read_id in event["reads"]:
            event_signatures.append(signatures[read_id])

        start_values = []
        end_values = []

        for signature in event_signatures:
            start_values.append(signature["overlap_start"])
            end_values.append(signature["overlap_end"])

        rows.append({
            "INVERSION": event["inversion"],
            "EVENT_ID": event_id,
            "EVENT_TYPE": "LOCAL_OVERLAP",
            "N_READS": len(event["reads"]),
            "READ_IDS": ";".join(event["reads"]),
            "IS_RECURRENT": len(event["reads"]) >= 2,
            "CHR": event["chr"],
            "ANCHOR_TYPE": event["anchor_type"],
            "MEDIAN_ANCHOR": median(event["anchor_values"]),
            "MEDIAN_START": median(start_values),
            "MEDIAN_END": median(end_values),
            "N_SA": 0,
            "MEDIAN_SA_START": None,
            "MEDIAN_SA_END": None
        })
        
    for event_id, event in gap_events.items():
        event_signatures = []
        
        for read_id in event["reads"]:
            event_signatures.append(signatures[read_id])

        start_values = []
        end_values = []
        sa_starts=[]
        sa_ends=[]
        for signature in event_signatures:
            start_values.append(signature["region_start"])
            end_values.append(signature["region_end"])
            for sa in signature["sa_support"]:
                sa_starts.append(sa["sa_start"])
                sa_ends.append(sa["sa_end"])
        
        rows.append({
            "INVERSION": event["inversion"],
            "EVENT_ID": event_id,
            "EVENT_TYPE": "GAP",
            "N_READS": len(event["reads"]),
            "READ_IDS": ";".join(event["reads"]),
            "IS_RECURRENT": len(event["reads"]) >= 2,
            "CHR": event["chr"],
            "ANCHOR_TYPE": event["anchor_type"],
            "MEDIAN_ANCHOR": median(event["anchor_values"]),
            "MEDIAN_START": median(start_values),
            "MEDIAN_END": median(end_values),
            "N_SA": len(sa_starts),
            "MEDIAN_SA_START": median(sa_starts),
            "MEDIAN_SA_END": median(sa_ends)
            })
    
    return pd.DataFrame(rows)

if __name__ == "__main__":

    RESULTS_DIR = Path("results/invdup_analysis")
    for bam_dir in RESULTS_DIR.iterdir():
        if bam_dir.is_dir():
            bam_name = bam_dir.name
            process_candidate_coordinates(
                bam_dir / f"{bam_name}_classification.tsv",
                bam_dir / f"{bam_name}_read_segments.tsv",
                bam_dir / f"{bam_name}_pairs.tsv"
            )


"""
OVERLAP:
vale la osa es que la 0231 te SUPER ESTABLE START LOCAL_OVERLAP
→ mateix chr
→ coordenada recurrent (start o end) dins ±COORDINATE_MARGIN
for cada signature:
    mirar events compatibles
    si n'hi ha:
        afegir read
    si no:
        crear event

GAP:
agrupar per reads de REGIONS i dins aquests mirar si tenen el mateix SA
mateixa inversió
mateix chr dels REGION
mateix strand dels REGION
flank_1 semblant ±COORDINATE_MARGIN
flank_2 semblant ±COORDINATE_MARGIN
i almenys un SA apunta a una zona semblant

requisit mida:
event_size = (
    sum(current_event["overlap_sizes"])
    / len(current_event["overlap_sizes"])
)

similar_size = abs(overlap_size - event_size) <= OVERLAP_SIZE_MARGIN
"""