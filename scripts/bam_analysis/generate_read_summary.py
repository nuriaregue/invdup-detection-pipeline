"""
1. Llegir allcoords
2. Llegir blast_pairs.tsv
3. Obrir cada BAM
4. Classificar reads amb utils_bam.process_reads()
5. Anotar SA amb find_supplimentary.search_supplementary()
6. Anotar self-BLAST amb process_overlapping.summarize_read_pair_overlaps()
7. Escriure read_summary.tsv
"""
from pathlib import Path
import csv

from scripts.bam_analysis import utils_bam as ub
from scripts.bam_analysis import find_supplementary as fs
from scripts.bam_analysis import process_overlapping as po

from scripts.config import (
    BAM_DIR,
    ALLCOORDS_PATH,
    MARGIN,
    BLAST_PAIRS_TSV,
    SUMMARY_TSV
)


POSS_PATTERNS = {
    "ONLY_SUPPL",
    "SPLIT_INV",
    "POSS_INV_DUP",
    "SPLIT"
}


def select_interesting_reads(reads_by_inv: dict) -> dict:
    """
    Keep only reads that are useful for the read-level summary. Simple reads are already ruled out.
    """
    interesting_reads = {}

    for read_id, read in reads_by_inv.items():
        if read["pattern"] in POSS_PATTERNS:
            interesting_reads[read_id] = read

    return interesting_reads


def decide_read(read: dict, overlap_summary: dict) -> str:
    """
    TODO: millors criteris
    Decide what to do with one read based on:
        - regional pattern
        - SA chromosome information
        - self-BLAST overlap ambiguity
    It used to obtain the reads to prioritize
    """
    decision = "LOW_PRIORITY"

    pattern = read.get("pattern", "")
    has_sa = read.get("has_sa_tag", False)
    sa_status = read.get("sa_recovery_status", "NO_SA")
    overlaps_both_blocks = overlap_summary.get("overlaps_both_blocks", False)

    if overlaps_both_blocks:
        decision = "REPEAT_PAIR_AMBIGUOUS"

    elif pattern in POSS_PATTERNS:
        decision = "KEEP_FOR_BREAKPOINT_ANALYSIS"

    elif pattern == "ONLY_SUPPL" and sa_status == "SA_OTHER_CHROMOSOMES":
        decision = "ONLY_SUPPL_WITH_OTHER_CHROMOSOME_SA"

    elif pattern == "ONLY_SUPPL" and sa_status == "SA_RECOVERED":
        decision = "KEEP_FOR_BREAKPOINT_ANALYSIS"

    elif pattern == "ONLY_SUPPL" and has_sa:
        decision = "ONLY_SUPPL_WITH_SAME_CHROMOSOME_SA"

    elif pattern == "ONLY_SUPPL":
        decision = "ONLY_SUPPL_NO_SA"

    return decision


def read_to_summary_row(bam_name: str,inv: str,read_id: str,read: dict,overlap_summary: dict,decision: str) -> dict:
    """
    Convert one read dictionary into one TSV row.
    """
    alignments=read.get("alignments",[])
    low_mapq,min_mapq=ub.classify_alignments_by_mapq_min(alignments) #we classify each alignment according to its mapq

    return {
        "BAM": bam_name,
        "INVERSION": inv,
        "READ_ID": read_id,
        "N_TOTAL_SEGMENTS": len(alignments),
        "N_SEGMENTS_GOOD_MAPQ":len(min_mapq),
        "N_SEGMENTS_LOW_MAPQ":len(low_mapq),
        "FLAGS_REGION": ",".join(map(str, read.get("flags", []))),
        "PATTERN_REGION": read.get("pattern", ""),
        "HAS_SA_TAG": read.get("has_sa_tag", False),
        "N_SA_TAGS": read.get("n_sa_tags", 0),
        "SA_CHROMOSOMES": ";".join(read.get("sa_chromosomes", [])),
        "SA_RECOVERY_STATUS": read.get("sa_recovery_status", "NA"),
        "N_RECOVERED_SA_ALIGNMENTS": read.get("n_recovered_sa_alignments", 0),
        "N_UNRECOVERED_SA_TAGS": read.get("n_unrecovered_sa_tags", 0),
        "SOME_SELF_BLAST_OVERLAP": overlap_summary.get("some_overlap", False),
        "OVERLAPS_BOTH_BLOCKS": overlap_summary.get("overlaps_both_blocks", False),
        "PAIR_IDS_BOTH_BLOCKS": ";".join(overlap_summary.get("pair_ids_both_blocks", [])),
        "DECISION": decision,
    }


def write_tsv(output_path: Path, rows: list) -> None:
    """
    Write all read summary rows to a TSV file.
    """
    fieldnames = [
        "BAM",
        "INVERSION",
        "READ_ID",
        "N_TOTAL_SEGMENTS",
        "N_SEGMENTS_GOOD_MAPQ",
        "N_SEGMENTS_LOW_MAPQ",
        "FLAGS_REGION",
        "PATTERN_REGION",
        "HAS_SA_TAG",
        "N_SA_TAGS",
        "SA_CHROMOSOMES",
        "SA_RECOVERY_STATUS",
        "N_RECOVERED_SA_ALIGNMENTS",
        "N_UNRECOVERED_SA_TAGS",
        "SOME_SELF_BLAST_OVERLAP",
        "OVERLAPS_BOTH_BLOCKS",
        "PAIR_IDS_BOTH_BLOCKS",
        "DECISION",
    ]

    with output_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def process_one_bam(bam_path: Path,all_inv: dict,repeated_regions: dict) -> list:
    """
    Process one BAM and return read_summary rows.
    """
    rows = []
    inv = ub.get_inv_from_bam_name(bam_path)

    if inv in all_inv:

        classified_reads = ub.process_reads(
            bam_path=bam_path,
            allcoords_inv={inv: all_inv[inv]},
            margin=MARGIN,
        )

        interesting_reads = select_interesting_reads(classified_reads[inv])

        if len(interesting_reads) == 0:
            print(f"No interesting reads found for {bam_path.name}")

        inv_chr = all_inv[inv][0]

        fs.search_supplementary(
            classified_reads=interesting_reads,
            bam_path=bam_path,
            inv_chr=inv_chr,
        )

        regions_in_inv = repeated_regions.get(inv, {})

        for read_id, read in interesting_reads.items():
            overlap_summary = po.summarize_read_pair_overlaps(read=read,regions_in_inv=regions_in_inv,min_overlap_bp=100,min_overlap_fraction=0.20,)
            decision = decide_read(read, overlap_summary)

            row = read_to_summary_row(
                bam_name=bam_path.name,
                inv=inv,
                read_id=read_id,
                read=read,
                overlap_summary=overlap_summary,
                decision=decision,
            )

            rows.append(row)

    else:
        print(f"Skipping {bam_path.name}: inversion not found in allcoords")

    return rows


if __name__ == "__main__":
    """
    Generate read_summary.tsv from:
        - BAM alignments
        - SA tag annotation
        - self-BLAST overlap annotation
    """
    

    all_inv = ub.proces_allcoords(ALLCOORDS_PATH, MARGIN)
    repeated_regions = po.generate_repeted_regions_dict(BLAST_PAIRS_TSV)

    rows = []

    bam_paths = sorted(BAM_DIR.glob("*.bam"))

    for bam_path in bam_paths:
        bam_rows = process_one_bam(
            bam_path=bam_path,
            all_inv=all_inv,
            repeated_regions=repeated_regions,
        )

        for row in bam_rows:
            rows.append(row)

    write_tsv(SUMMARY_TSV, rows)

    print(f"Written: {SUMMARY_TSV}")
    print(f"Number of rows: {len(rows)}")

