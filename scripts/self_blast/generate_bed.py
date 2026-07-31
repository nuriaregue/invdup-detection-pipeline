import csv
from pathlib import Path

from scripts.config import (
    BLAST_PAIRS_TSV,
    BED_DIR,
)


def needs_update(input_path: Path, output_path: Path) -> bool:
    """
    Return True if the output file does not exist or is older than the input file.
    """

    if not output_path.exists():
        return True

    return input_path.stat().st_mtime > output_path.stat().st_mtime


def read_pairs_by_inversion(pairs_tsv: Path) -> dict:
    """
    Read blast_pairs.tsv and group pair rows by inversion.
    """

    pairs_by_inv = {}

    with pairs_tsv.open("r", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:
            inv = row["INVERSION"]

            if inv not in pairs_by_inv:
                pairs_by_inv[inv] = []

            pairs_by_inv[inv].append(row)

    return pairs_by_inv


def pair_to_bed_rows(pair: dict) -> list:
    """
    Convert one deduplicated repeat pair into two BED rows:
    one row for block A and one row for block B.
    """

    chrom = pair["CHROMOSOME"]

    block_a_start = int(pair["BLOCK_A_START"])
    block_a_end = int(pair["BLOCK_A_END"])
    block_b_start = int(pair["BLOCK_B_START"])
    block_b_end = int(pair["BLOCK_B_END"])

    pair_id = pair["PAIR_ID"]

    bed_rows = [
        {
            "chrom": chrom,
            "chromStart": block_a_start - 1,
            "chromEnd": block_a_end,
            "name": f"{pair_id}_A",
        },
        {
            "chrom": chrom,
            "chromStart": block_b_start - 1,
            "chromEnd": block_b_end,
            "name": f"{pair_id}_B",
        },
    ]

    return bed_rows


def create_bed_rows(pairs: list) -> list:
    """
    Convert all pairs from one inversion into BED rows.
    """

    bed_rows = []

    for pair in pairs:
        bed_rows.extend(pair_to_bed_rows(pair))

    return bed_rows


def write_bed_file(bed_rows: list, output_path: Path) -> None:
    """
    Write BED rows to a BED file.
    """

    fieldnames = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
    ]

    with output_path.open("w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writerows(bed_rows)


def generate_bed(pairs_tsv: Path, bed_dir: Path) -> None:
    """
    Generate one BED file per inversion from blast_pairs.tsv.
    BED files are regenerated only if they are missing or older than blast_pairs.tsv.
    """

    bed_dir.mkdir(parents=True, exist_ok=True)

    pairs_by_inv = read_pairs_by_inversion(pairs_tsv)

    for inv, pairs in pairs_by_inv.items():
        bed_path = bed_dir / f"{inv}_region.bed"

        if needs_update(pairs_tsv, bed_path):
            print(f"Creating/updating {bed_path}...")
            bed_rows = create_bed_rows(pairs)
            write_bed_file(bed_rows, bed_path)
        else:
            print(f"{bed_path} already exists and is up to date.")


if __name__ == "__main__":
    generate_bed(
        pairs_tsv=BLAST_PAIRS_TSV,
        bed_dir=BED_DIR,
    )