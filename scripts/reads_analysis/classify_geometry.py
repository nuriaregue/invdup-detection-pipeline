
"""
1. Mirar la geometria del read --> Com estan organitzats els segments dins d'aquest read?
2. Mirar si aquest patró es repeteix en la resta de reads? Els reads que tenen una geometria semblant també connecten les mateixes regions del genoma?


"""
from pathlib import Path
from scripts.config import (
    MIN_GEOMETRY_MAPQ,
    MIN_SA_RECOVERY_MAPQ,
    MIN_GAP_COVERAGE_FRACTION
)

import pandas as pd
from pandas import DataFrame

def process_segments_tsv(segments_tsv:Path,output_path: Path,pairs_path:Path)->DataFrame:
    if not segments_tsv.exists():
        raise FileNotFoundError(segments_tsv)
    else:
        summary_rows=[]
        all_pair_rows=[]
        df = pd.read_csv(segments_tsv,delimiter="\t")
        #Each group contains all alignment segments associated with one read.
        total_reads = df.groupby(["BAM", "READ_ID"])
        for key, read_df in total_reads:
            bam = key[0]
            read_id = key[1]
            #We keep all alignment segments to preserve the complete read-level metadata and alignment counts.
            all_alignments = read_df.sort_values(by="ALIGNMENT_INDEX")
            #We calculate geometric relationships only from alignments that pass the minimum mapping-quality threshold.
            geometry_alignments = read_df[read_df["MAPQ"] >= MIN_GEOMETRY_MAPQ].copy()
            #We keep reliable supplementary alignments as supporting evidence even when they do not pass the stricter geometry MAPQ threshold.
            suplementary_alignments=read_df[read_df["SOURCE"].isin(["RECOVERED_SA", "SA_TAG_ONLY"]) & (read_df["MAPQ"] >= MIN_SA_RECOVERY_MAPQ)]
            #We order the selected segments according to their corrected position along the original read.
            geometry_alignments = geometry_alignments.sort_values(by=["CORR_QUERY_START","CORR_QUERY_END","ALIGNMENT_INDEX"])
            suplementary_alignments = suplementary_alignments.sort_values(by=["CORR_QUERY_START","CORR_QUERY_END","ALIGNMENT_INDEX"])
            row,pair_rows = process_one_read(bam,read_id,all_alignments,geometry_alignments,suplementary_alignments)
            summary_rows.append(row)  #Each read produces one summary row.
            all_pair_rows.extend(pair_rows) #We add each row individually to the list with extends

    summary_df = pd.DataFrame(summary_rows)
    summary_df.columns=summary_df.columns.str.upper()
    summary_df.to_csv(str(output_path),sep="\t",index=False)

    pairs_df = pd.DataFrame(all_pair_rows)
    pairs_df.columns=pairs_df.columns.str.upper()
    pairs_df.to_csv(str(pairs_path),sep="\t",index=False)
    return summary_df

def calculate_merged_coverage(intervals: list):
    """
    Calculate the total number of query bases covered by one or more
    intervals without counting overlapping bases more than once.
    """
    if not intervals:
        return 0
    intervals = sorted(intervals)

    current_start = intervals[0][0]
    current_end = intervals[0][1]
    total_coverage = 0
    for start, end in intervals[1:]:
        # Overlapping intervals are merged into the same covered region.
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total_coverage += current_end - current_start
            current_start = start
            current_end = end
    total_coverage += current_end - current_start
    return total_coverage

def has_sa_in_region_gap(geom_align: DataFrame,suplementary_alignments: DataFrame):
    first_region,second_region,supporting_sa,gap_type = find_sa_supported_region_gap(
        geom_align,
        suplementary_alignments
    )

    return gap_type

def find_sa_supported_region_gap(geom_align,suplementary_alignments):
    regions = geom_align[geom_align["SOURCE"] == "REGION"].sort_values(by=["CORR_QUERY_START","CORR_QUERY_END"]).to_dict("records")

    for i in range(len(regions)-1):
        first_region = regions[i]
        second_region = regions[i+1]

        gap_start = first_region["CORR_QUERY_END"]
        gap_end = second_region["CORR_QUERY_START"]

        same_chromosome = first_region["REF_CHR"] == second_region["REF_CHR"]
        same_strand = first_region["STRAND"] == second_region["STRAND"]

        if gap_end > gap_start and same_chromosome and same_strand:
            gap_length = gap_end-gap_start

            supporting_intervals = []
            supporting_chromosomes = []
            supporting_sa = []
            has_opposite_sa = False

            for sa in suplementary_alignments.to_dict("records"):
                overlap_start = max(gap_start,sa["CORR_QUERY_START"])
                overlap_end = min(gap_end,sa["CORR_QUERY_END"])

                if overlap_end > overlap_start:
                    supporting_intervals.append((overlap_start,overlap_end))
                    supporting_chromosomes.append(sa["REF_CHR"])
                    supporting_sa.append(sa)

                    if sa["STRAND"] != first_region["STRAND"]:
                        has_opposite_sa = True

            covered_bp = calculate_merged_coverage(supporting_intervals)
            coverage_fraction = covered_bp/gap_length

            if has_opposite_sa and coverage_fraction >= MIN_GAP_COVERAGE_FRACTION:
                if any(chromosome != first_region["REF_CHR"] for chromosome in supporting_chromosomes):
                    gap_type = "INTERCHROMOSOMAL"
                else:
                    gap_type = "INTRACHROMOSOMAL"

                return first_region,second_region,supporting_sa,gap_type

    return None,None,[],""

def build_initial_row(bam: str,read_id: str,all_align: DataFrame,geom_align: DataFrame,) -> dict:
    #General read metadata and alignment counts
    n_segments_total = len(all_align)  # Total count includes every alignment segment of the read
    n_geometry_segments = len(geom_align) # Only segments that pass the geometry MAPQ threshold
    n_recovered_segments = int((geom_align["SOURCE"] == "RECOVERED_SA").sum()) #Count recovered alignments
    return {
        "bam": bam,
        "inversion": all_align["INVERSION"].iloc[0],
        "read_id": read_id,
        "pattern_region": all_align["PATTERN_REGION"].iloc[0],
        "sa_recovery_status": all_align["SA_RECOVERY_STATUS"].iloc[0],
        "n_segments_total": n_segments_total,
        "n_geometry_segments": n_geometry_segments,
        "n_low_mapq_segments": n_segments_total - n_geometry_segments,
        "n_region_segments": int((geom_align["SOURCE"] == "REGION").sum()),
        "n_recovered_segments": n_recovered_segments,
        "has_recovered_sa": n_recovered_segments > 0,
        #Patterns describe only the segments selected for geometry
        "source_pattern": ";".join(geom_align["SOURCE"].astype(str)), 
        "flag_pattern": ";".join(geom_align["FLAG"].apply(lambda flag: "NA" if pd.isna(flag) else str(int(flag)))),
        "ref_intervals": "",
        "strand_pattern": "",
        "query_gaps_bp": [],
        "query_overlaps_bp": [],
        "max_query_gap_bp": None,
        "max_query_overlap_bp": None,
        "max_ref_overlap_bp": None,
        "min_mapq": None,
        "has_sa_in_region_gap": False,
        "sa_region_gap_type": ""
    }

def add_geometry(row_dict: dict,geom_align: DataFrame, pair_rows: list,suplementary_alignments:DataFrame) -> dict:
    """
    Add segment-level information and summarize the pair-level metrics
    calculated for the read.
    """
    #No reliable alignments available
    if geom_align.empty: #If no alignment passes the mapq filter, we return base dict
        return row_dict
    read_dict = geom_align.to_dict("records") #The list is converted into a list of dictionaries to process each alignment
    
    #We first store the reference interval, strand and MAPQ of the available geometry segments. These values do not
    # require a comparison between two segments
    row_dict["ref_intervals"] = get_ref_intervals(read_dict)
    row_dict["strand_pattern"] = get_full_pattern(read_dict)
    row_dict["min_mapq"] = get_min_mapq(read_dict)

    #Gap and overlap geometry requires at least two alignments
    if not pair_rows:
        return row_dict

    #Query gaps and overlaps have already been calculated for each and consecutive pair, so we reuse them for the read-level summary.
    query_gaps = []
    for pair in pair_rows:
        if pair["IS_CONSECUTIVE"]:
            query_gaps.append(pair["QUERY_GAP_BP"])
    query_overlaps = []
    for pair in pair_rows:
        if pair["IS_CONSECUTIVE"]:
            query_overlaps.append(pair["QUERY_OVERLAP_BP"])
    #Reference overlap is calculated separately because it compares all possible segment pairs, not only consecutive pairs
    ref_overlaps = process_ref_overlap(read_dict) 
    #The full lists are preserved, while the maximum values provide a simple read-level summary
    gap_type=has_sa_in_region_gap(geom_align,suplementary_alignments)
    sa_has_gap=False
    if gap_type !="":
        sa_has_gap=True
    #Update basic dict
    row_dict.update({
        "query_gaps_bp": query_gaps,
        "query_overlaps_bp": query_overlaps,
        "max_query_gap_bp": max(query_gaps, default=0),
        "max_query_overlap_bp": max(query_overlaps, default=0),
        "max_ref_overlap_bp": max(ref_overlaps, default=0),
        "has_sa_in_region_gap": sa_has_gap, #We use reliable supplementary alignments to check whether they provide opposite-strand support inside a gap between regional segments.
        "sa_region_gap_type":gap_type
    })

    return row_dict

def process_segment_pairs(bam: str,read_id: str,inversion: str,geom_align: DataFrame) -> list:
    """
    Generate one row for every possible geometry segment pair.
    """
    pair_rows = []
    #We order the segments according to their corrected position along the original read.
    geom_align = geom_align.sort_values(by=["CORR_QUERY_START","CORR_QUERY_END","ALIGNMENT_INDEX"])
    alignments = geom_align.to_dict("records")
    pair_number = 0

    for i in range(len(alignments)):
        for j in range(i + 1, len(alignments)):
            first_align = alignments[i]
            second_align = alignments[j]
            same_chromosome = first_align["REF_CHR"] == second_align["REF_CHR"]
            opposite_strands = first_align["STRAND"] != second_align["STRAND"]

            query_gap,query_overlap =  calculate_query_relationship(first_align,second_align)

            ref_results = calculate_reference_relationship(first_align,second_align,same_chromosome)
            ref_overlap_bp = ref_results[0]
            reference_overlap_fraction = ref_results[1]
            reference_length_ratio = ref_results[2]

            pair_rows.append({
                "BAM": bam,
                "INVERSION": inversion,
                "READ_ID": read_id,
                "PAIR_INDEX": pair_number,
                "ALIGNMENT_INDEX_A": first_align["ALIGNMENT_INDEX"],
                "ALIGNMENT_INDEX_B": second_align["ALIGNMENT_INDEX"],
                "IS_CONSECUTIVE": j == i + 1,
                "SAME_CHROMOSOME": same_chromosome,
                "OPPOSITE_STRANDS": opposite_strands,
                "QUERY_GAP_BP": query_gap,
                "QUERY_OVERLAP_BP": query_overlap,
                "REFERENCE_OVERLAP_BP": ref_overlap_bp,
                "REFERENCE_OVERLAP_FRACTION": reference_overlap_fraction,
                "REFERENCE_LENGTH_RATIO": reference_length_ratio
            })
            pair_number += 1
            
    return pair_rows

def calculate_reference_relationship(current: dict, next_segment: dict, same_chromosome: bool):
    """
    Calculate the reference overlap and compare the segment lengths.
    """
    ref_len_a = current["REF_END"] - current["REF_START"] 
    ref_len_b = next_segment["REF_END"] - next_segment["REF_START"]

    short_len = min(ref_len_a, ref_len_b)
    long_len = max(ref_len_a, ref_len_b)

    if same_chromosome:
        overlap_start = max(current["REF_START"], next_segment["REF_START"])
        overlap_end = min(current["REF_END"], next_segment["REF_END"])
        ref_overlap = max(0, overlap_end - overlap_start)
    else:
        ref_overlap = None

    if ref_overlap is not None and short_len > 0:
        overlap_fraction_short = ref_overlap / short_len
    else:
        overlap_fraction_short = None

    if long_len > 0:
        segment_length_ratio = short_len / long_len
    else:
        segment_length_ratio = None

    return ref_overlap, overlap_fraction_short, segment_length_ratio

def process_one_read(bam:str,read_id:str,all_align:DataFrame,geom_align:DataFrame,suplementary_alignments:DataFrame):
    """
    Generate the read-level summary and the pair-level row for one read
    """
    row_dict = build_initial_row(bam,read_id,all_align,geom_align)
    # We first calculate the metrics of each consecutive segment pair.
    pair_rows = process_segment_pairs(bam,read_id,all_align["INVERSION"].iloc[0],geom_align)
    #The read-level summary reuses the pair-level query metrics.
    row_dict = add_geometry(row_dict=row_dict,geom_align=geom_align,pair_rows=pair_rows,suplementary_alignments=suplementary_alignments)
    return row_dict, pair_rows

def calculate_query_relationship(current: dict,next_segment: dict):
    """
    Calculate the query gap or overlap between two consecutive segments.
    """
    query_diff = next_segment["CORR_QUERY_START"] - current["CORR_QUERY_END"]
    
    #A >0 difference represents an uncovered gap between the segments,while a <0 difference means that their query intervals overlap.
    if query_diff < 0:
        query_gap = 0
        query_overlap = abs(query_diff)
    else:
        query_gap = query_diff
        query_overlap = 0

    return query_gap, query_overlap

def process_ref_overlap(read:list):
    """
    Calculate reference overlaps between alignment segments located on the
    same chromosome and in opposite orientations.
    """
    ref_overlaps = []
    for i in range(len(read)):
        current_alignment = read[i]
        for j in range(i+1,len(read)):
            next_alignment = read[j]
            #Reference overlap is calculated only for segments that map to the same chromosome and have opposite orientations.
            are_same_chr=current_alignment["REF_CHR"]==next_alignment["REF_CHR"]
            are_opp_strand=current_alignment["STRAND"] != next_alignment["STRAND"]
            if are_opp_strand and are_same_chr: 
                """
                To calculate the shared reference interval, we take the latest start coordinate because the overlap cannot begin
                before both segments have started.
                We then take the earliest end coordinate because the overlap finishes as soon as one of the two segments ends.
                If the earliest end is greater than the latest start,the difference between them represents the overlap length
                """
                start_overlap=max(current_alignment["REF_START"],next_alignment["REF_START"])
                end_overlap=min(current_alignment["REF_END"],next_alignment["REF_END"])
                overlap=end_overlap-start_overlap
                if overlap>0:
                    ref_overlaps.append(overlap)
                else:
                    ref_overlaps.append(0)
    return ref_overlaps


def get_full_pattern(read_dict:list):
    """
    Build the strand pattern of the read by concatenating the orientation
    of its geometry alignment segments.
    """
    strand_pattern = ""
    for alignment in read_dict:
        strand_pattern += alignment["STRAND"]
    return strand_pattern

def get_min_mapq(read_dict:list):
    """
    Return the lowest mapping-quality value among the alignment segments
    used for geometry.
    """
    mapq_values=[]
    for alignment in read_dict:
        mapq_values.append(alignment["MAPQ"])
    return min(mapq_values,default=0)

def get_ref_intervals(read:list):
    """
    Store the reference chromosome and coordinate interval of each geometry
    alignment segment in a single string.
    """
    ref_intervals = []
    for alignment in read:
        interval = f"{alignment['REF_CHR']}:{alignment['REF_START']}-{alignment['REF_END']}"
        ref_intervals.append(interval)
    return ";".join(ref_intervals)

def summarize_geometry(summary_df: DataFrame) -> DataFrame:
    """
    Generate a BAM-level summary from the read geometry table.
    """
    if summary_df.empty:
        return pd.DataFrame()

    summary = summary_df.groupby("BAM").agg(
        N_READS=("READ_ID", "nunique"),
        MEDIAN_QUERY_GAP_BP=("MAX_QUERY_GAP_BP", "median"),
        MEAN_QUERY_GAP_BP=("MAX_QUERY_GAP_BP", "mean"),
        MEDIAN_QUERY_OVERLAP_BP=("MAX_QUERY_OVERLAP_BP", "median"),
        MEAN_QUERY_OVERLAP_BP=("MAX_QUERY_OVERLAP_BP", "mean"),
        MEDIAN_REFERENCE_OVERLAP_BP=("MAX_REF_OVERLAP_BP", "median"),
        MEAN_REFERENCE_OVERLAP_BP=("MAX_REF_OVERLAP_BP", "mean")
    ).reset_index()

    return summary
