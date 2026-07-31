"""
Generate a first read-level INV_DUP candidate classification from the geometry
and consecutive segment-pair tables.
"""

from pathlib import Path

import pandas as pd
from pandas import DataFrame


# Provisional thresholds. We keep them here so they can be adjusted after
# reviewing the first candidate distributions.
MAX_QUERY_GAP_BP = 200
MAX_QUERY_OVERLAP_BP = 200
MIN_MIRROR_REF_OVERLAP_FRACTION = 0.80
MIN_MIRROR_LENGTH_RATIO = 0.80


def has_mirror_like_pair(read_pairs: DataFrame) -> bool:
    """
    Check whether one consecutive pair has a strongly mirror-like geometry.
    """
    if read_pairs.empty:
        return False
    overlap_column = "REFERENCE_OVERLAP_FRACTION"
    length_column = "REFERENCE_LENGTH_RATIO"

    if overlap_column not in read_pairs.columns or length_column not in read_pairs.columns:
        return False

    for pair in read_pairs.itertuples():
        same_chromosome = bool(pair.SAME_CHROMOSOME)
        opposite_strands = bool(pair.OPPOSITE_STRANDS)
        overlap_fraction = pair.REFERENCE_OVERLAP_FRACTION
        length_ratio = pair.REFERENCE_LENGTH_RATIO

        values_available = (
            pd.notna(overlap_fraction)
            and pd.notna(length_ratio)
        )
        is_consecutive=pair.IS_CONSECUTIVE
        if values_available:
            if (same_chromosome and opposite_strands and is_consecutive and
            overlap_fraction >= MIN_MIRROR_REF_OVERLAP_FRACTION and length_ratio >= MIN_MIRROR_LENGTH_RATIO):
                return True

    return False

def classify_one_read(row: pd.Series,read_pairs: DataFrame) -> tuple:
    """
    Assign a simple provisional classification and the main reason.
    """
    n_segments = int(row.get("N_GEOMETRY_SEGMENTS",0))
    strand_pattern = str(row.get("STRAND_PATTERN",""))
    max_query_gap = float(row.get("MAX_QUERY_GAP_BP",0))
    max_query_overlap = float(row.get("MAX_QUERY_OVERLAP_BP",0))
    ref_intervals = str(row.get("REF_INTERVALS",""))

    #We require at least one repeated chromosome in reads with three or more segments, although coordinate compatibility is checked later.
    has_repeated_chr = has_repeated_chromosome(ref_intervals)
    #An alignment recovered from or represented by the SA tag may explain a query gap that is not covered by the REGION alignments.
    has_sa_in_region_gap = row.get("HAS_SA_IN_REGION_GAP",False)
    sa_region_gap_type = str( row.get("SA_REGION_GAP_TYPE", ""))

    #Classification:
    if n_segments < 2: # Fewer than two reliable segments are not enough to compare geometry. 
        return "NOT_ENOUGH_INFORMATION","INSUFFICIENT_GEOMETRY"
    elif strand_pattern == "": # An empty strand pattern indicates that the geometry could not be evaluated.
        return "AMBIGUOUS","NO_STRAND_PATTERN"
    elif high_percentage_overlap(max_query_overlap):  #A large query overlap is considered incompatible because the two alignments may represent the same part of the read.
        return "NOT_INV_DUP","LARGE_QUERY_OVERLAP"
    
    # Main geometry with only one orientation does not support an inversion.
    # However, an opposite-strand supplementary alignment inside the gap may complete the inverted pattern.
    elif len(set(strand_pattern)) == 1: #TODO:revisar aquest cas
        if has_sa_in_region_gap: 
            return "POSSIBLE_INV_DUP","REGION_GAP_WITH_INVERTED_SA_SUPPORT"
        elif max_query_gap > MAX_QUERY_GAP_BP: 
            return "AMBIGUOUS","SAME_STRAND_WITH_UNEXPLAINED_QUERY_GAP"
        else:
            return "NOT_INV_DUP","SAME_STRAND_ONLY"
        
    elif has_mirror_like_pair(read_pairs): # Strongly symmetric opposite-strand mappings are preserved as ambiguous because they may represent redundant or uncertain alignments.
        return "AMBIGUOUS","MIRROR_LIKE_GEOMETRY" #We discard the whole read
    
    elif n_segments >= 3:
        # Reads with three or more segments require at least one repeated chromosome to preserve the genomic context of the event.
        if not has_repeated_chr:
            return "NOT_INV_DUP","NO_REPEATED_CHROMOSOME"
        # Opposite-strand segments from the same chromosome with reference overlap support a local intrachromosomal INV_DUP.
        elif has_inverted_ref_overlap(read_pairs):
            return "POSSIBLE_INV_DUP","COMPATIBLE_INTRACHROMOSOMAL_OVERLAP"
        # An inverted SA alignment from another chromosome explains the original query gap between the regional flanks.
        elif sa_region_gap_type == "INTERCHROMOSOMAL":
            return "POSSIBLE_INV_DUP","INTERCHROMOSOMAL_GAP_WITH_INVERTED_SA"
        # An inverted SA alignment from a distant region of the same chromosome explains the original query gap between the regional flanks.
        elif sa_region_gap_type == "INTRACHROMOSOMAL":
            return "POSSIBLE_INV_DUP","DISTANT_INTRACHROMOSOMAL_GAP_WITH_INVERTED_SA"
        # A large query gap is rejected when no supplementary alignment explains the missing sequence.
        elif max_query_gap > MAX_QUERY_GAP_BP:
            return "NOT_INV_DUP","LARGE_UNEXPLAINED_QUERY_GAP"
        # Mixed-strand geometries without reference overlap or supplementary gap support are preserved for later analysis.
        else:
            return "AMBIGUOUS","INVERTED_PATTERN_WITHOUT_ENOUGH_SUPPORT"
        
    else: #We have also reads with exactly two reliable geometry segments, that could be INV_DUP
        #A large unexplained query gap makes the two-segment geometry insufficiently continuous.
        if max_query_gap > MAX_QUERY_GAP_BP and not has_sa_in_region_gap:
            return "NOT_INV_DUP","LARGE_UNEXPLAINED_QUERY_GAP"
        #Two opposite-strand segments overlapping on the same reference chromosome provide partial evidence of an inverted duplication.
        elif has_inverted_ref_overlap(read_pairs):
            return "POSSIBLE_INV_DUP","COMPATIBLE_TWO_SEGMENT_PATTERN"
        #Opposite-strand segments without reference overlap do not match the current expected two-segment INV_DUP geometry.
        else:
            return "NOT_INV_DUP","INVERTED_TWO_SEGMENTS_NON_OVERLAP"


def has_inverted_ref_overlap(read_pairs: DataFrame) -> bool:
    """
    Check whether an opposite-strand pair overlaps in reference coordinates.
    """
    inverted_ref_overlaps=read_pairs[
        (read_pairs["SAME_CHROMOSOME"]==True) & (read_pairs["OPPOSITE_STRANDS"]==True) & (read_pairs["REFERENCE_OVERLAP_BP"] > 0)]
    return not inverted_ref_overlaps.empty

def high_percentage_overlap(query_ov: float) -> bool:
    """
    Detect long read overlaps, they can be found in inversions, but not in INV_DUP.
    """
    return query_ov >= MAX_QUERY_OVERLAP_BP

def has_repeated_chromosome(ref_intervals: str) -> bool:
    """
    Check whether at least one chromosome appears more than once.
    Examples:
        chr9:100-200;chr1:500-700;chr9:300-400 -> True
        chr12:100-200;chr10:500-700;chr16:300-400 -> False
    """
    chromosomes = []

    # Each reference interval has the format chromosome:start-end.We keep only the chromosome name from every alignment.
    intervals = ref_intervals.split(";")
    for interval in intervals:
        chromosome = interval.split(":")[0]
        chromosomes.append(chromosome)

    # If the number of chromosomes is greater than the number of unique chromosomes, at least one chromosome appears more than once.
    unique_chromosomes = set(chromosomes)
    if len(chromosomes) > len(unique_chromosomes):
        return True
    else:
        return False
def process_geometry_tsv(geometry_path: Path,pairs_path: Path,output_path: Path) -> DataFrame:
    """
    Read the geometry tables, classify every read and write candidate_reads.tsv.
    """
    if not geometry_path.exists():
        raise FileNotFoundError(geometry_path)

    geometry_df = pd.read_csv(geometry_path, sep="\t")

    if pairs_path.exists() and pairs_path.stat().st_size > 0:
        pairs_df = pd.read_csv(pairs_path, sep="\t")
    else:
        pairs_df = pd.DataFrame()

    candidate_rows = []

    for _, row in geometry_df.iterrows():
        bam = row["BAM"]
        read_id = row["READ_ID"]

        if pairs_df.empty:
            read_pairs = pd.DataFrame()
        else:
            read_pairs = pairs_df[(pairs_df["BAM"] == bam) & (pairs_df["READ_ID"] == read_id)]

        candidate_type, reason = classify_one_read(row, read_pairs)

        candidate_row = row.to_dict()
        candidate_row["CANDIDATE_TYPE"] = candidate_type
        candidate_row["CLASSIFICATION_REASON"] = reason
        candidate_row["HAS_MIRROR_LIKE_PAIR"] = has_mirror_like_pair(read_pairs)

        candidate_rows.append(candidate_row)

    candidate_df = pd.DataFrame(candidate_rows)
    candidate_df.to_csv(output_path, sep="\t", index=False)

    return candidate_df


def summarize_candidates(candidate_df: DataFrame) -> DataFrame:
    """
    Count reads by candidate type and classification reason.
    """
    if candidate_df.empty:
        return pd.DataFrame()
    
    return candidate_df.groupby(["BAM", "CANDIDATE_TYPE", "CLASSIFICATION_REASON"])["READ_ID"].nunique().reset_index(name="N_READS")
