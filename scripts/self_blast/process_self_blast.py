from pathlib import Path
import csv
from scripts.self_blast import deduplicate_self_blast as dd
from scripts.self_blast import run_self_blast as rsb

from scripts.config import (
    BLAST_RESULTS_DIR,
    OUTPUT_TSV_DIR,
    ALLCOORDS_PATH,
    BLAST_HITS_TSV,
    BLAST_PAIRS_TSV,
    BLAST_SUMMARY_TSV,
    MARGIN,
    MIN_IDENTITY,
    MIN_LENGTH,
)

"""
Process self-BLAST hit tables for inversion regions.

This script:
1. Reads BLAST outfmt 6 tables from BLASTS/BLAST_RESULTS.
2. Filters strong non-self hits.
3. Summarizes hits per inversion.
4. Writes blast_summary.tsv and blast_hits.tsv.
"""
def is_self_hit(qstart: int, qend: int, sstart: int, send: int) -> bool:
    """
    Check if the BLAST hit maps the same query interval to the same subject interval.
    """

    return qstart == sstart and qend == send


def get_blast_orientation(sstart: int, send: int) -> str:
    """
    Return the orientation of the BLAST hit based on subject coordinates.
    """

    if sstart < send:
        return "direct"

    return "inverted"


def is_not_hidden_file(filename: str) -> bool:
    """
    Ignore hidden files when scanning folders.
    """

    return not filename.startswith(".")


def find_hit_tables(path: Path) -> list:
    """
    Find BLAST hit table result files in the specified folder.
    """
    hit_tables = []

    for file_path in path.iterdir():
        if file_path.is_file() and is_not_hidden_file(file_path.name) and file_path.suffix == ".txt":
            hit_tables.append(file_path)
    return hit_tables


def get_inv_id_from_blast_file(file_path: Path) -> str:
    """
    Extract inversion ID from BLAST result filename.

    Example:
    ./BLASTS/BLAST_RESULTS/HsInv0036_region.txt -> HsInv0036
    """
    inv_id = file_path.stem #stem is filename without extension
    inv_id = inv_id.replace("_region", "")
    return inv_id


def filter_blast_results(file_path: Path,min_identity: float = 90.0,min_length: int = 1000) -> list:
    """
    Filter BLAST outfmt 6 hits using identity, length and self-hit filters.
    """

    if not file_path.exists():
        raise FileNotFoundError(f"{file_path} was not found!")

    filtered_hits = []

    with file_path.open("r") as blast_result:
        #BLAST hit tables are 1-based
        for line in blast_result:
            columns = line.strip().split("\t")

            pident = float(columns[2])
            length = int(columns[3])
            qstart = int(columns[6])
            qend = int(columns[7])
            sstart = int(columns[8])
            send = int(columns[9])

            self_hit = is_self_hit(qstart, qend, sstart, send)
            orientation = get_blast_orientation(sstart, send)

            if pident >= min_identity and length >= min_length and not self_hit:
                hit = {
                    "pident": pident,
                    "length": length,
                    "qstart": qstart,
                    "qend": qend,
                    "sstart": sstart,
                    "send": send,
                    "orientation": orientation,
                }

                filtered_hits.append(hit)

    return filtered_hits



def print_blast_summary(inv_id: str, summary: dict) -> None:
    """
    Print a short summary of the filtered BLAST hits for one inversion.
    """

    print(f"INVERSION: {inv_id}")

    if summary["n_hits_filtered"] > 0:
        print(f"\t n_hits_filtered: {summary['n_hits_filtered']}")
        print(f"\t n_direct_hits: {summary['n_direct_hits']}")
        print(f"\t n_inverted_hits: {summary['n_inverted_hits']}")
        print(f"\t max_hit_length: {summary['max_hit_length']}")
        print(f"\t max_hit_identity: {summary['max_hit_identity']}")
    else:
        print("\t No hits meet the established filters.")


def process_blast_result(file_path: Path) -> dict:
    """
    Method to process the self-BLAST result for one inversion.
    """
    inv_id = get_inv_id_from_blast_file(file_path)
    hits = filter_blast_results(file_path,MIN_IDENTITY,MIN_LENGTH)
    summary = dd.summarize_blast_hits(hits)
    print_blast_summary(inv_id, summary)
    return {
        inv_id: {
            "source_file": file_path,
            "pair_counter": 1,
            "params": {
                "min_identity": MIN_IDENTITY,
                "min_length": MIN_LENGTH,
                "margin": MARGIN
            },
            "summary": summary,
            "hits": hits
        }
    }


def write_blast_summary_tsv(results_by_inv: dict, output_path: Path):
    """
    Write one summary row per inversion to a TSV file.
    """

    fieldnames = [
        "INVERSION",
        "N_HITS_FILTERED",
        "N_UNIQUE_PAIRS",
        "N_DIRECT_HITS",
        "N_UNIQUE_DIRECT_PAIRS",
        "N_INVERTED_HITS",
        "N_UNIQUE_INVERTED_PAIRS",
        "MAX_HIT_LENGTH",
        "MAX_HIT_IDENTITY",
        "MIN_IDENTITY",
        "MIN_LENGTH",
        "MARGIN",
    ]

    rows = []

    for inv_id, data in results_by_inv.items():
        summary = data["summary"]
        summary_pairs=data["summary_pairs"]
        params = data["params"]

        row = {
            "INVERSION": inv_id,
            "N_HITS_FILTERED": summary["n_hits_filtered"],
            "N_UNIQUE_PAIRS":summary_pairs["n_unique_pairs"],
            "N_DIRECT_HITS": summary["n_direct_hits"],
            "N_UNIQUE_DIRECT_PAIRS":summary_pairs["n_unique_direct_pairs"],
            "N_INVERTED_HITS": summary["n_inverted_hits"],
            "N_UNIQUE_INVERTED_PAIRS":summary_pairs["n_unique_inverted_pairs"],
            "MAX_HIT_LENGTH": summary["max_hit_length"],
            "MAX_HIT_IDENTITY": summary["max_hit_identity"],
            "MIN_IDENTITY": params["min_identity"],
            "MIN_LENGTH": params["min_length"],
            "MARGIN": params["margin"]
        }

        rows.append(row)

    with output_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)


def write_blast_hits_tsv(results_by_inv: dict, output_path: Path):
    """
    Write one row per filtered BLAST hit to a TSV file.
    """

    fieldnames = [
        "INVERSION",
        "PIDENT",
        "LENGTH",
        "QSTART",
        "QEND",
        "SSTART",
        "SEND",
        "ORIENTATION"]

    rows = []

    for inv_id, data in results_by_inv.items():
        for hit in data["hits"]:
            row = {
                "INVERSION": inv_id,
                "PIDENT": hit["pident"],
                "LENGTH": hit["length"],
                "QSTART": hit["qstart"],
                "QEND": hit["qend"],
                "SSTART": hit["sstart"],
                "SEND": hit["send"],
                "ORIENTATION": hit["orientation"],
            }

            rows.append(row)

    with output_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)
    
