"""
Run the INV_DUP analysis pipeline for selected BAM files.

For each BAM, the pipeline:
1. Selects candidate reads.
2. Recovers reliable supplementary alignments.
3. Generates the segment-level TSV.
4. Generates the read geometry TSV.
"""

from pathlib import Path
import pandas as pd
import sys
from scripts.config import (
    BAM_DIR,
    ALLCOORDS_PATH,
    MARGIN,
    PIPELINE_INV_DUP_RESULTS_DIR,
    MIN_READS_INV_DUP,
    FASTA_DIR,
    BLAST_RESULTS_DIR,
    OUTPUT_TSV_DIR,
    REFERENCE_PATH,
    BLAST_HITS_TSV,
    BLAST_PAIRS_TSV,
    BLAST_SUMMARY_TSV,
    BED_DIR,
    MIN_IDENTITY,
    MIN_LENGTH
)

from scripts.bam_analysis import utils_bam as ub
from scripts.bam_analysis import generate_read_segments as grs
from scripts.reads_analysis import classify_geometry as cg
from scripts.reads_analysis import classify_inv as ci
from scripts.reads_analysis import classify_coords as cc
from scripts.self_blast import run_self_blast as rsb
from scripts.self_blast import process_self_blast as psb
from scripts.self_blast import deduplicate_self_blast as dd 
from scripts.self_blast import generate_bed as gb
from scripts.reads_analysis import compare_selfblast as cs

"""
    BAM_DIR / "HsInv0170_HG00599.bam",
    BAM_DIR / "HsInv0170_HG00639.bam",
    BAM_DIR / "HsInv0186_HG03563.bam",
    BAM_DIR / "HsInv0186_NA20355.bam",
    BAM_DIR / "HsInv1146_HG02258.bam",
    BAM_DIR / "HsInv1146_HG02280.bam",
    BAM_DIR / "HsInv0001_Std.bam",
    BAM_DIR / "HsInv1153_HG02841.bam",
"""

def process_self_blast(inv_id: str, bam_name: str):
    """
    Method to generate and process the reference self-BLAST results.
    """
    print("\n" + "*" * 60)
    print(f"CREATING SELF-BLAST FOR {inv_id} WITH {MIN_IDENTITY}% min_identity AND {MIN_LENGTH} bp as min_length")
    print("*" * 60)

    #The self-BLAST FASTA and raw BLAST result directories are shared between BAMs.
    FASTA_DIR.mkdir(parents=True, exist_ok=True)
    BLAST_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    #Each BAM has its own folder for the processed self-BLAST TSVs.
    blast_output_dir = OUTPUT_TSV_DIR / bam_name
    blast_output_dir.mkdir(parents=True, exist_ok=True)

    blast_hits_path = blast_output_dir / f"{bam_name}_blast_hits.tsv"
    blast_pairs_path = blast_output_dir / f"{bam_name}_blast_pairs.tsv"
    blast_summary_path = blast_output_dir / f"{bam_name}_blast_summary.tsv"

    #We select only the inversion associated with the current BAM.
    allcoords_inv = rsb.read_inversions_from_allcoords(ALLCOORDS_PATH)
    inversions = {inv_id: allcoords_inv[inv_id]}

    rsb.process_index(REFERENCE_PATH)
    rsb.generate_fasta_blast(REFERENCE_PATH,inversions,MARGIN,FASTA_DIR)

    fasta_path = FASTA_DIR / f"{inv_id}_region.fa"
    hit_table_path = rsb.generate_blast_hit_table(fasta_path,BLAST_RESULTS_DIR)
    results_by_inv = psb.process_blast_result(hit_table_path)

    psb.write_blast_hits_tsv(results_by_inv,blast_hits_path)
    dd.write_blast_pairs_tsv(results_by_inv,blast_pairs_path,inversions)
    psb.write_blast_summary_tsv(results_by_inv,blast_summary_path)

    bed_output_dir = BED_DIR / bam_name
    bed_output_dir.mkdir(parents=True, exist_ok=True)
    gb.generate_bed(blast_pairs_path,bed_output_dir)

def process_bam(bam_path: Path,all_inv: dict) -> None:
    """
    Run the complete read-segment and geometry analysis for one BAM.
    """
    print("\n" + "*" * 85)
    print(f"Processing: {bam_path.name}")
    print("*" * 85)

    if not bam_path.exists():
        print(f"BAM not found: {bam_path}")
    else:
        segment_rows = grs.process_one_bam(bam_path=bam_path,all_inv=all_inv)
        if len(segment_rows) == 0:
            print("No candidate segments found")
        else:
            bam_name=bam_path.stem #With stem, extension (.bam) is removed
            bam_output_dir = PIPELINE_INV_DUP_RESULTS_DIR / bam_name
            bam_output_dir.mkdir(parents=True,exist_ok=True)

            inv_id = bam_name.split("_")[0]
            #process_self_blast(inv_id,bam_name)
            blast_pairs_path = OUTPUT_TSV_DIR / bam_name / f"{bam_name}_blast_pairs.tsv"

            print("\n" + "*" * 60)
            print(f"FILTERING INV_DUP FOR {bam_name}")
            print("*" * 60)

            #Each TSV stores the analysis at a different level:
            # - read_segments: individual alignments
            # - segments_geometry: one geometry summary per read
            # - segment_pairs: geometry and information between consecutive reliable alignments
            segments_path = bam_output_dir / f"{bam_name}_read_segments.tsv"
            geometry_path = bam_output_dir / f"{bam_name}_segments_geometry.tsv"
            pairs_path=bam_output_dir / f"{bam_name}_pairs.tsv"
            candidates_path=bam_output_dir / f"{bam_name}_classification.tsv"
            grs.write_tsv(output_path=segments_path,rows=segment_rows)

            print(f"Segment rows written: {len(segment_rows)}")
            print(f"Segments TSV: {segments_path}")

            geometry_df = cg.process_segments_tsv(segments_tsv=segments_path,output_path=geometry_path,pairs_path=pairs_path)

            candidate_df = ci.process_geometry_tsv(geometry_path=geometry_path,pairs_path=pairs_path,output_path=candidates_path)
            signatures, events_df = cc.process_candidate_coordinates(candidates_path,segments_path,pairs_path)
            bam_summary = summarize_bam(candidate_df,events_df)
            summary_df = pd.DataFrame([bam_summary])
            summary_path = bam_output_dir / f"{bam_name}_summary.tsv"
            summary_df.to_csv(summary_path,sep="\t",index=False)
            # We compare the main candidate event with self-BLAST pairs already generated for this BAM.
            if blast_pairs_path.exists():
                cs.filter_reference_repeats(summary_path,blast_pairs_path)
            else:
                print(f"Self-BLAST pairs not found: {blast_pairs_path}")

            print(f"Candidate rows written: {len(candidate_df)}")
            print(f"Candidates TSV: {candidates_path}")
            if candidate_df.empty:
                print("No geometry rows generated")
            else:
                geometry_summary = ci.summarize_candidates(candidate_df)
                print("\nClassification summary:")
                print(geometry_summary.to_string(index=False))

            if (not candidate_df.empty and "HAS_RECOVERED_SA" in candidate_df.columns):
                n_recovered_reads = int(geometry_df["HAS_RECOVERED_SA"].sum())
                print(f"Reads with recovered SA: {n_recovered_reads}")

            # The final summary includes the read classification, coordinate event and reference-repeat validation.
            important_fields = [
                "BAM", "INVERSION", "N_READS", "N_POSSIBLE_INV_DUP", "N_EVENTS", "N_RECURRENT_EVENTS",
                "INV_DUP_EVENT_TYPE", "INV_DUP_MAIN_EVENT_READS", "EVENT_START", "EVENT_END", "FRACTION_POSSIBLE_READS_IN_MAIN_EVENT",
                "IS_POSS_INV_DUP", "HAS_REFERENCE_REPEAT", "REFERENCE_PAIR_ID", "REFERENCE_ORIENTATION", "REFERENCE_FILTER_STATUS"
            ]
            final_summary_df = pd.read_csv(summary_path,delimiter="\t")

            print("\n" + "*" * 60)
            print(f"FINAL SUMMARY FOR {bam_name}")
            print("*" * 60)
            print(final_summary_df[important_fields].T.to_string(header=False))
            print("\nConsult the final TSV to see the supporting READ_IDs.")

def summarize_bam(candidate_df:pd.DataFrame, events_df:pd.DataFrame) -> dict:
    """
    Method to summarize read classifications and coordinate support for one BAM.
    """
    n_reads = candidate_df["READ_ID"].nunique()
    possible_df = candidate_df[candidate_df["CANDIDATE_TYPE"] == "POSSIBLE_INV_DUP"]
    n_possible = possible_df["READ_ID"].nunique()
    n_not_inv_dup = candidate_df[candidate_df["CANDIDATE_TYPE"] == "NOT_INV_DUP"]["READ_ID"].nunique()
    n_ambiguous = candidate_df[candidate_df["CANDIDATE_TYPE"] == "AMBIGUOUS"]["READ_ID"].nunique()
    possible_read_ids = possible_df["READ_ID"].tolist()
    n_not_enough = candidate_df[candidate_df["CANDIDATE_TYPE"] == "NOT_ENOUGH_INFORMATION"]["READ_ID"].nunique()

    n_events = len(events_df)
    if events_df.empty:
        n_recurrent_events = 0
    else:
        recurrent_events = events_df[events_df["IS_RECURRENT"] == True]
        n_recurrent_events = len(recurrent_events)

    main_event = None

    if not events_df.empty:
        main_event = events_df.loc[
            events_df["N_READS"].idxmax()
        ]

    summary = {
        "BAM": candidate_df["BAM"].iloc[0],
        "INVERSION": candidate_df["INVERSION"].iloc[0],
        "N_READS": n_reads,
        "N_POSSIBLE_INV_DUP": n_possible,
        "N_NOT_INV_DUP": n_not_inv_dup,
        "N_NOT_ENOUGH_INFORMATION": n_not_enough,
        "N_AMBIGUOUS": n_ambiguous,
        "POSSIBLE_READ_IDS": ";".join(possible_read_ids),
        "N_EVENTS": n_events,
        "N_RECURRENT_EVENTS": n_recurrent_events
    }

    if main_event is not None:
        summary["INV_DUP_EVENT_ID"] = main_event["EVENT_ID"]
        summary["INV_DUP_EVENT_TYPE"] = main_event["EVENT_TYPE"]
        summary["INV_DUP_MAIN_EVENT_READS"] = main_event["N_READS"]
        summary["INV_DUP_EVENT_READ_IDS"] = main_event["READ_IDS"]
        summary["INV_DUP_EVENT_CHR"] = main_event["CHR"]
        summary["EVENT_START"] = main_event["MEDIAN_START"]
        summary["EVENT_END"] = main_event["MEDIAN_END"]

        if main_event["EVENT_TYPE"] == "GAP":
            summary["N_SA"] = main_event["N_SA"]
            summary["SA_MEDIAN_START"] = main_event["MEDIAN_SA_START"]
            summary["SA_MEDIAN_END"] = main_event["MEDIAN_SA_END"]
        else:
            summary["N_SA"] = 0
            summary["SA_MEDIAN_START"] = None
            summary["SA_MEDIAN_END"] = None

        if n_possible > 0:
            summary["FRACTION_POSSIBLE_READS_IN_MAIN_EVENT"] = (
                main_event["N_READS"] / n_possible
            )
        else:
            summary["FRACTION_POSSIBLE_READS_IN_MAIN_EVENT"] = 0

        is_recurrent=main_event["IS_RECURRENT"]
        summary["HAS_RECURRENT_COORDINATES"] = is_recurrent
        if is_recurrent and n_possible >=MIN_READS_INV_DUP:
            summary["IS_POSS_INV_DUP"]="YES"
        elif n_possible<MIN_READS_INV_DUP and is_recurrent:
            summary["IS_POSS_INV_DUP"]="NOT_ENOUGH_READS"
        else:
            summary["IS_POSS_INV_DUP"]="NO"
    else:
        summary["INV_DUP_EVENT_ID"] = ""
        summary["INV_DUP_EVENT_TYPE"] = ""
        summary["INV_DUP_MAIN_EVENT_READS"] = 0
        summary["INV_DUP_EVENT_READ_IDS"] = ""
        summary["INV_DUP_EVENT_CHR"] = ""
        summary["EVENT_START"] = None
        summary["EVENT_END"] = None
        summary["N_SA"] = 0
        summary["SA_MEDIAN_START"] = None
        summary["SA_MEDIAN_END"] = None
        summary["FRACTION_POSSIBLE_READS_IN_MAIN_EVENT"] = 0
        summary["HAS_RECURRENT_COORDINATES"] = False
        summary["IS_POSS_INV_DUP"] = "NO"
    return summary


if __name__ == "__main__":
    PIPELINE_INV_DUP_RESULTS_DIR.mkdir(parents=True,exist_ok=True)

    if len(sys.argv) != 2:
        print("USAGE: python3 -m scripts.run_pipeline BAM_NAME.bam")
        sys.exit(1)

    bam_name = sys.argv[1]
    bam_path = BAM_DIR / bam_name

    all_inv = ub.proces_allcoords(ALLCOORDS_PATH,MARGIN)
    process_bam(bam_path=bam_path,all_inv=all_inv)

"""
TODO:
1. Setup del environment més automàtic + script per posar nom del bam
2. Comentaris
3. Taula resum 
4. Revisar tot el tema dels comentaris
5. Distribució del 75% com fa el ricardo (si es canvia dir-ho)
6. Imatges!!!
8. Explicar que és cada tsv i la classificacio final
"""