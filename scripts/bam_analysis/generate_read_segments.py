from pathlib import Path
import csv

from scripts.bam_analysis import utils_bam as ub
from scripts.bam_analysis import find_supplementary as fs

from scripts.config import (
    BAM_DIR,
    ALLCOORDS_PATH,
    MARGIN,
    SEGMENTS_TSV
)


CANDIDATE_PATTERNS = {
    "POSS_INV_DUP",
    "SPLIT", 
    "SPLIT_INV",
    "ONLY_SUPPL"
}

def get_segment_sort_key(row: dict):
    """
    Return the corrected query coordinates used to sort segments.
    """
    return row["CORR_QUERY_START"], row["CORR_QUERY_END"]

def get_hard_clip_start(cigar: str, strand: str) -> int:
    """
    Method that returns the number of hard-clipped bases located before the
    alignment in the original read orientation.
    For forward alignments, this corresponds to the beginning of the CIGAR.
    For reverse alignments, this corresponds to the end of the CIGAR.
    Examples:
        cigar = 100H500M, strand = "+" -> 100
        cigar = 500M100H, strand = "-" -> 100
        cigar = 100H500M, strand = "-" -> 0
    """
    if cigar is None:
        return 0
    if strand == "+":
        number = ""
        for char in cigar:
            if char.isdigit():
                number += char
            else:
                if char == "H" and number != "":
                    return int(number)
                else:
                    return 0
    elif strand == "-":
        #For reverse alignments, the clipping before the segment in the
        #original read appears at the end of the CIGAR.
        number = ""
        if not cigar.endswith("H"):
            return 0
        index = len(cigar) - 2 # Skip the final H and start at the preceding digit
        while index >= 0 and cigar[index].isdigit():
            number = cigar[index] + number
            index -= 1     
        if number != "":
            return int(number)
        else:
            return 0
    else:
        raise ValueError(
            f"Invalid strand '{strand}'. Expected '+' or '-'."
        )

    return 0

def select_candidate_reads(reads_by_inv: dict) -> dict:
    """
    Method to keep reads that are useful for segment-level analysis.
    """
    candidate_reads = {}

    for read_id, read in reads_by_inv.items():
        if read["pattern"] in CANDIDATE_PATTERNS:
            candidate_reads[read_id] = read

    return candidate_reads


def get_strand(alignment) -> str:
    """
    Method that returns strand from BAM flag.
    """
    if alignment.is_reverse:
        return "-"
    else:
        return "+"

def process_query_coordinates(query_start:int,query_end:int, query_length: int,cigar:str,strand:str):
    """
    Method that returns the correct query alignment coordinates
    """
    # Some supplementary alignments start with hard clipping in the CIGAR.
    # In these cases, pysam query_alignment_start can be 0 because the hard-clipped
    # bases are not stored in that BAM record. To recover the position in the
    # original read, we add hard-clipped length according to its strand orientation 
    # (soft-clipped are already stored in the BAM record)
    coord={}

    hard_clip_start = get_hard_clip_start(cigar,strand)
    if strand == "+":
        coord["start"] = query_start + hard_clip_start
        coord["end"] = query_end + hard_clip_start
    elif strand == "-":
        coord["start"] = query_length - query_end + hard_clip_start
        coord["end"] = query_length - query_start + hard_clip_start
    else:
        raise ValueError(f"Invalid strand '{strand}'. Expected '+' or '-'.")
    return coord

def alignment_to_segment_row(bam_name: str,inv: str,read_id: str,read: dict,alignment,source: str) -> dict:
    """
    Method to convert one alignment into one row of read_segments.tsv.
    Source describes where the segment was obtained from, not its BAM alignment type. Primary/supplementary/secondary status is stored in separate columns.
    """
    if alignment.query_length is None:
        raise ValueError(f"Query length is unavailable for read {read_id}, flag {alignment.flag}, CIGAR {alignment.cigarstring}")
    strand=get_strand(alignment)
    #Convert pysam query coordinates to the original read orientation.
    #Reverse alignments are reflected and hard clipping is restored.
    query_coord = process_query_coordinates(alignment.query_alignment_start,alignment.query_alignment_end,
                                            alignment.query_length,alignment.cigarstring,strand)
    return {
        "BAM": bam_name,
        "INVERSION": inv,
        "READ_ID": read_id,
        "SOURCE": source, 
        "FLAG": alignment.flag,
        "STRAND": strand,
        "REF_CHR": alignment.reference_name,
        "REF_START": alignment.reference_start,
        "REF_END": alignment.reference_end,
        "QUERY_START": alignment.query_alignment_start,
        "CORR_QUERY_START":query_coord["start"],
        "QUERY_END": alignment.query_alignment_end,
        "CORR_QUERY_END":query_coord["end"],
        "QUERY_ALIGNED_BP":  alignment.query_alignment_length,
        "MAPQ": alignment.mapping_quality,
        "CIGAR": alignment.cigarstring,
        "IS_PRIMARY": not alignment.is_secondary and not alignment.is_supplementary,
        "IS_SUPPLEMENTARY": alignment.is_supplementary,
        "IS_SECONDARY": alignment.is_secondary,
        "PATTERN_REGION": read.get("pattern", ""),
        "FLAGS_REGION": ";".join(map(str, read.get("flags", []))),
        "SA_RECOVERY_STATUS": read.get("sa_recovery_status", "NA"),
    }

def sa_tag_to_segment_row(bam_name: str, inv: str, read_id: str, read: dict, sa_tag: dict) -> dict:
    """
    Convert an unrecovered SA tag into one row of segments.tsv.

    Important:
        - QUERY_START and QUERY_END are left empty because the complete alignment
        record is not available in the current subset BAM. For comparisons, we
        use CORR_QUERY_START and CORR_QUERY_END inferred from the SA CIGAR.
        - SA tags do not contain the complete BAM flag, so I cannot know if the
        alignment is primary, supplementary or secondary. Empty values mean
        unknown, not False.
    """
    cigar = sa_tag["cigar"]
    strand = sa_tag["strand"]

    reference_start = sa_tag["start"] - 1 # SA POS is 1-based, while pysam and the output TSV use 0-based coordinates.
    reference_length = ub.get_reference_length_from_cigar(cigar)
    reference_end = reference_start + reference_length

    query_coordinates = ub.get_sa_query_coordinates(cigar, strand)

    return {
        "BAM": bam_name,
        "INVERSION": inv,
        "READ_ID": read_id,
        "SOURCE": "SA_TAG_ONLY",
        "FLAG": "",
        "STRAND": strand,
        "REF_CHR": sa_tag["chr"],
        "REF_START": reference_start,
        "REF_END": reference_end,
        "QUERY_START": "",
        "CORR_QUERY_START": query_coordinates["start"],
        "QUERY_END": "",
        "CORR_QUERY_END": query_coordinates["end"],
        "QUERY_ALIGNED_BP": ub.get_query_alignment_length_from_cigar(cigar),
        "MAPQ": sa_tag["mapq"],
        "CIGAR": cigar,
        "IS_PRIMARY": "",
        "IS_SUPPLEMENTARY": "",
        "IS_SECONDARY": "",
        "PATTERN_REGION": read.get("pattern", ""),
        "FLAGS_REGION": ";".join(map(str, read.get("flags", []))),
        "SA_RECOVERY_STATUS": read.get("sa_recovery_status", "NA")
    }

def read_to_segment_rows(bam_name: str,inv: str,read_id: str,read: dict,include_recovered_sa: bool = False) -> list:
    """
    Method to convert all alignments of one read into segment rows.
    Alignments in our exploring region are always included.
    SA-recovered alignments are included only when include_recovered_sa is True.
    """
    rows = []

    alignments = read.get("alignments", [])

    for alignment in alignments:
        row = alignment_to_segment_row(bam_name,inv,read_id,read,alignment,"REGION")
        rows.append(row)

    if include_recovered_sa:
        for alignment in read.get("recovered_sa_alignments", []):
            row = alignment_to_segment_row(bam_name,inv,read_id,read,alignment,"RECOVERED_SA")
            rows.append(row)

    for sa_tag in read.get("unrecovered_sa_tags", []):
        row = sa_tag_to_segment_row(bam_name, inv, read_id, read, sa_tag)
        rows.append(row)

    #Order segments according to their position in the original read.
    rows.sort(key=get_segment_sort_key)
    #Add an alignment_index, once the alignments are sorted by their corrected coordinates
    alignment_index = 0
    for row in rows:
        row["ALIGNMENT_INDEX"] = alignment_index
        alignment_index += 1
    return rows


def write_tsv(output_path: Path, rows: list) -> None:
    """
    Method to write read segment rows to TSV.
    """
    fieldnames = [
        "BAM",
        "INVERSION",
        "READ_ID",
        "ALIGNMENT_INDEX",
        "SOURCE",
        "FLAG",
        "STRAND",
        "REF_CHR",
        "REF_START",
        "REF_END",
        "QUERY_START",
        "CORR_QUERY_START",
        "QUERY_END",
        "CORR_QUERY_END",
        "QUERY_ALIGNED_BP",
        "MAPQ",
        "CIGAR",
        "IS_PRIMARY",
        "IS_SUPPLEMENTARY",
        "IS_SECONDARY",
        "PATTERN_REGION",
        "FLAGS_REGION",
        "SA_RECOVERY_STATUS"
    ]

    with output_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def process_one_bam(bam_path: Path, all_inv: dict) -> list:
    """
    Method to generate segment rows for one BAM.
    """
    rows = []

    inv = ub.get_inv_from_bam_name(bam_path)

    if inv not in all_inv:
        print(f"Skipping {bam_path.name}: inversion not found in allcoords")
        return rows

    classified_reads = ub.process_reads(bam_path=bam_path,allcoords_inv={inv: all_inv[inv]},margin=MARGIN)

    candidate_reads = select_candidate_reads(classified_reads[inv])

    inv_chr = all_inv[inv][0]

    fs.search_supplementary(classified_reads=candidate_reads,bam_path=bam_path,inv_chr=inv_chr)

    for read_id, read in candidate_reads.items():
        segment_rows = read_to_segment_rows(bam_name=bam_path.name,inv=inv,read_id=read_id,read=read,include_recovered_sa=True)

        for row in segment_rows:
            rows.append(row)

    return rows


if __name__ == "__main__":
    all_inv = ub.proces_allcoords(ALLCOORDS_PATH, MARGIN)

    rows = []

    bam_paths = sorted(BAM_DIR.glob("*.bam"))

    for bam_path in bam_paths:
        bam_rows = process_one_bam(bam_path=bam_path,all_inv=all_inv)

        for row in bam_rows:
            rows.append(row)

    write_tsv(SEGMENTS_TSV, rows)

    print(f"Written: {SEGMENTS_TSV}")
    print(f"Number of segment rows: {len(rows)}")

