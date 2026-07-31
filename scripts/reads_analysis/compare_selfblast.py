from pathlib import Path
import pandas as pd

from scripts.config import COORDINATE_MARGIN


def coordinates_are_close(position_1: int, position_2: int) -> bool:
    """
    Method to check if two genomic positions are within the coordinate margin.
    """
    return abs(position_1 - position_2) <= COORDINATE_MARGIN


def find_reference_repeat(summary_df: pd.DataFrame, blast_pairs_path: Path) -> dict:
    """
    Method to check if the main INV_DUP event could be explained by a repeat already present in the reference.
    """
    result = {
        "HAS_REFERENCE_REPEAT": False,
        "REFERENCE_PAIR_ID": "",
        "REFERENCE_ORIENTATION": "",
        "REFERENCE_FILTER_STATUS": "KEEP"
    }

    blast_pairs_df = pd.read_csv(blast_pairs_path,delimiter="\t")

    # The comparison is only possible when a main candidate event and self-BLAST pairs are available.
    if not blast_pairs_df.empty and summary_df["INV_DUP_EVENT_ID"].iloc[0] != "":
        event_chr = summary_df["INV_DUP_EVENT_CHR"].iloc[0]
        event_start = summary_df["EVENT_START"].iloc[0]
        event_end = summary_df["EVENT_END"].iloc[0]

        for pair in blast_pairs_df.itertuples():
            same_chromosome = pair.CHROMOSOME == event_chr

            # We focus on inverted reference repeats because they can reproduce an INV_DUP-like alignment pattern.
            inverted_repeat = pair.ORIENTATION == "inverted"
            block_a_start = pair.BLOCK_A_START - 1
            block_a_end = pair.BLOCK_A_END
            block_b_start = pair.BLOCK_B_START - 1
            block_b_end = pair.BLOCK_B_END

            block_a_start_close = coordinates_are_close(block_a_start,event_start)
            block_a_end_close = coordinates_are_close(block_a_end,event_end)

            block_b_start_close = coordinates_are_close(block_b_start,event_start)
            block_b_end_close = coordinates_are_close(block_b_end,event_end)

            # The event overlaps coordinates belonging to one of the inverted repeat blocks in the reference.
            block_a_match = block_a_start_close or block_a_end_close
            block_b_match = block_b_start_close or block_b_end_close

            if same_chromosome and inverted_repeat and (block_a_match or block_b_match):
                result["HAS_REFERENCE_REPEAT"] = True
                result["REFERENCE_PAIR_ID"] = pair.PAIR_ID
                result["REFERENCE_ORIENTATION"] = pair.ORIENTATION
                result["REFERENCE_FILTER_STATUS"] = "REVIEW_REFERENCE_REPEAT"

    return result


def filter_reference_repeats(summary_path: Path, blast_pairs_path: Path):
    """
    Method to compare the main candidate event with repeats detected in the reference.
    """
    summary_df = pd.read_csv(summary_path,delimiter="\t")

    reference_repeat = find_reference_repeat(summary_df,blast_pairs_path)

    # We keep the read-based classification and add the reference-repeat evidence as a final validation layer.
    summary_df["HAS_REFERENCE_REPEAT"] = reference_repeat["HAS_REFERENCE_REPEAT"]
    summary_df["REFERENCE_PAIR_ID"] = reference_repeat["REFERENCE_PAIR_ID"]
    summary_df["REFERENCE_ORIENTATION"] = reference_repeat["REFERENCE_ORIENTATION"]
    summary_df["REFERENCE_FILTER_STATUS"] = reference_repeat["REFERENCE_FILTER_STATUS"]

    summary_df.to_csv(summary_path,sep="\t",index=False)

    return summary_df