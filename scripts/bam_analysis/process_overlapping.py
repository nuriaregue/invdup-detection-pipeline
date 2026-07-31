
from pathlib import Path
from scripts.bam_analysis import utils_bam as ub
from scripts.config import (
    BAM_DIR,
    ALLCOORDS_PATH,
    MARGIN,
    BLAST_PAIRS_TSV
)


def generate_repeted_regions_dict(path:Path):
    dict_regions={}
    with path.open() as pairs_tsv:
        header = pairs_tsv.readline().strip().split("\t")
        for line in pairs_tsv:
            if line.strip():
                columns=line.strip("\n").split("\t")
                inv=columns[1]
                if dict_regions.get(inv) is None:
                    dict_regions[inv] = {}
                #Convert 1-based hit table coordinates to 0-based coordinates compatible with pysam.
                dict_regions[inv][columns[0]] = {"block_a":{"start":int(columns[4])-1, "end":int(columns[5])},
                                                  "block_b":{"start":int(columns[6])-1, "end":int(columns[7])}, 
                                                  "orientation":columns[8], "chr":columns[9]}

                #we create a entry in the dictionary for each pair --> orientation?
                #list of regions? []
    
    return dict_regions


def get_interval(coord: dict):
    """
    Return start and end from a coordinate dictionary.
    """
    return coord["start"], coord["end"]

def interval_overlap(align_coords:dict, block_coords:dict):
    """
    Calculate the overlap between one alignment and one repeat block.

    Returns overlap_bp and overlap_perc.
    """
    start_align,end_align=get_interval(align_coords)
    start_block,end_block=get_interval(block_coords)
    start_overlap=max(start_block,start_align) 
    end_overlap=min(end_block,end_align)
    overlap_bp=end_overlap-start_overlap
    if overlap_bp<=0:
        return 0,0
    else:
        alignment_lenght=end_align-start_align
        overlap_perc=overlap_bp/alignment_lenght
        return overlap_bp,overlap_perc


def get_alignment_overlaps(alignment_coord: dict, regions_in_inv: dict, min_overlap_bp: int = 100, min_overlap_fraction: float = 0.20) -> list:
    """
    Find all repeat blocks overlapped by one alignment.

    For each overlap, it stores pair_id, block A/B, orientation,
    overlap_bp and overlap_fraction.
    """
    overlaped_regions=[]
    for pair_id,pair in regions_in_inv.items():
        block_a= pair.get("block_a",{})
        block_b= pair.get("block_b",{})
        overlap_bp_a,overlap_frac_a=interval_overlap(alignment_coord,block_a)
        overlap_bp_b,overlap_frac_b=interval_overlap(alignment_coord,block_b)
        if overlap_bp_a>=min_overlap_bp and overlap_frac_a>=min_overlap_fraction:
            #Block A overlaps with the alignment
            reg= {"pair_id":pair_id,"block":"A","orientation":pair.get("orientation"),
                                "overlap_bp":overlap_bp_a,"overlap_fraction":overlap_frac_a}
            overlaped_regions.append(reg)
        if overlap_bp_b>=min_overlap_bp and overlap_frac_b>=min_overlap_fraction:
            #Block B overlaps with the alignment
            reg= {"pair_id":pair_id,"block":"B","orientation":pair.get("orientation"),
                                "overlap_bp":overlap_bp_b,"overlap_fraction":overlap_frac_b}
            overlaped_regions.append(reg)
    return overlaped_regions

def init_read_summary(read: dict) -> dict:
    """
    Create an empty summary for one read.
    """
    alignments = read.get("alignments", [])
    read_id = None
    if len(alignments) > 0:
        read_id = alignments[0].query_name

    return {
        "read_id": read_id,"n_alignments": len(alignments),
        "pairs": {},
        "some_overlap":False,
        "overlaps_both_blocks": False,
        "pair_ids_both_blocks": []
    }

def add_overlap_to_summary(summary: dict, overlap: dict, alignment, alignment_index: int) -> None:
    """
    Add one alignment-block overlap to the read summary.
    """
    pair_id = overlap["pair_id"]
    block = overlap["block"]
    if pair_id not in summary["pairs"]:
        summary["pairs"][pair_id] = {
            "orientation": overlap.get("orientation"),
            "blocks_touched": set(),
            "touches_both_blocks": False,
            "overlaps": []
        }
    summary["pairs"][pair_id]["blocks_touched"].add(block)
    summary["pairs"][pair_id]["overlaps"].append({
        "alignment_index": alignment_index,
        "flag": alignment.flag,
        "ref_start": alignment.reference_start,
        "ref_end": alignment.reference_end,
        "query_start": alignment.query_alignment_start,
        "query_end": alignment.query_alignment_end,
        "block": block,
        "overlap_bp": overlap["overlap_bp"],
        "overlap_fraction": overlap["overlap_fraction"]
    })
def pair_touches_both_blocks(pair_info: dict) -> bool:
    """
    Return True if one pair has both block A and block B touched.
    """
    blocks = pair_info.get("blocks_touched", set())
    return "A" in blocks and "B" in blocks

def finalize_read_summary(summary: dict) -> dict:
    """
    Mark pairs where the read touches both blocks.
    """
    for pair_id, pair_info in summary["pairs"].items():
        if pair_touches_both_blocks(pair_info):
            pair_info["touches_both_blocks"] = True
            summary["pair_ids_both_blocks"].append(pair_id)
    if len(summary["pair_ids_both_blocks"]) > 0:
        summary["overlaps_both_blocks"] = True
    if len(summary["pairs"])>0:
        summary["some_overlap"]=True

    return summary

def summarize_read_pair_overlaps(read: dict,regions_in_inv: dict,min_overlap_bp: int = 100,min_overlap_fraction: float = 0.20) -> dict:
    """
    Summarize repeat-block overlaps for all alignments of one read.
    """
    summary = init_read_summary(read)
    alignments = read.get("alignments", [])
    alignment_index=0
    for alignment in alignments:
        alignment_coord = {
            "start": alignment.reference_start,
            "end": alignment.reference_end
        }
        overlaps = get_alignment_overlaps(alignment_coord,regions_in_inv,min_overlap_bp,min_overlap_fraction)
        for overlap in overlaps:
            add_overlap_to_summary(summary=summary,overlap=overlap,alignment=alignment,alignment_index=alignment_index)
        alignment_index+=1
    summary = finalize_read_summary(summary)

    return summary


def print_overlaping(classified_read: dict,summary:dict) -> str:
    both_overlapping=summary["overlaps_both_blocks"]
    some_overlap=summary["some_overlap"]
    alignments=classified_read.get("alignments",[])
    if len(alignments) > 0:
        read_id = alignments[0].query_name
    print(f"\t {read_id}|{classified_read['n_alignments']} alignments |" 
          f"flags = {classified_read['flags']}|pattern = {classified_read['pattern']}|"
          f"some_overlap = {some_overlap} | ignore = {both_overlapping}"
          )

if __name__ == "__main__":
    
    all_inv = ub.proces_allcoords(ALLCOORDS_PATH, MARGIN)
    repeated_regions = generate_repeted_regions_dict(BLAST_PAIRS_TSV)

    for bam_path in sorted(BAM_DIR.glob("*.bam")): #search for .bam files in the directory
        inv = ub.get_inv_from_bam_name(bam_path)
        print(f"\nProcessing BAM: {bam_path.name}")
        print(f"Inversion from BAM name: {inv}")

        classified_reads = ub.process_reads(bam_path=bam_path,allcoords_inv={inv: all_inv[inv]},margin=MARGIN)

        regions_in_inv = repeated_regions.get(inv, {})
        complex_reads = ub.get_reads_by_pattern(classified_reads[inv],"POSS_INV_DUP")
        print(f"Number of POSS_INV_DUP reads: {len(complex_reads)}")
        for read in complex_reads.values():
            summary = summarize_read_pair_overlaps(read, regions_in_inv,100,0.4)
            print_overlaping(read, summary)

   