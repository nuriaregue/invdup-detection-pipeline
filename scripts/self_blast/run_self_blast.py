from pathlib import Path

import pysam
import subprocess

"""
Generate self-BLAST input FASTA files and BLAST hit tables for inversion regions.

This script:
1. Reads inversion coordinates from BAM/allcoords.txt.
2. Extracts one FASTA region per inversion from the reference genome.
3. Runs BLAST of each extracted region against itself.
"""

from scripts.config import (
    ALLCOORDS_PATH,
    REFERENCE_PATH,
    FASTA_DIR,
    BLAST_RESULTS_DIR,
    MARGIN,
)


def read_inversions_from_allcoords(path:Path) -> dict:
    """
    Read inversion coordinates from allcoords.txt.
    """

    if not path.exists():
        raise FileNotFoundError(f"{path} was not found!")

    all_inv = {}

    with path.open("r") as allcoords:
        header = allcoords.readline().strip().split("\t")

        for line in allcoords:
            if line.strip():
                columns = line.strip().split("\t")
                # Match chromosome naming used in the reference FASTA.
                chromosome = "chr" + columns[1]
                all_inv[columns[0]] = {
                    "chr": chromosome,
                    "start": int(columns[2]),
                    "end": int(columns[5]),
                }

    return all_inv


def process_index(reference_path: Path):
    """
    Create FASTA index if it does not exist.
    """

    if not reference_path.exists():
        raise FileNotFoundError(f"{reference_path} was not found!")
    index_path=reference_path.with_name(reference_path.name+".fai")
    if not index_path.exists():
        print(f"FASTA index not found for {reference_path}! Creating index...")
        pysam.faidx(str(reference_path))
    


def generate_fasta_blast(reference_path: Path,invs: dict,margin: int,fasta_dir:Path):
    """
    Extract one FASTA region per inversion with a margin around it.
    """

    fasta_dir.mkdir(parents=True,exist_ok=True)

    with pysam.FastaFile(str(reference_path)) as ref_file:
        for inv_id, values in invs.items():

            start = values["start"] - margin
            end = values["end"] + margin

            if start < 1:
                start = 1

            # fetch() uses 0-based start and end-exclusive coordinates.
            fasta_str = ref_file.fetch(values["chr"], start - 1, end) #Is the same as samtools faidx Std_ref2.fa chr:start-end
            file_path = fasta_dir / f"{inv_id}_region.fa"

            try:
                with file_path.open("x") as output_fasta:
                    # FASTA header stores the genomic region used for the BLAST.
                    output_fasta.write(f">{values['chr']}:{start}-{end}\n")
                    output_fasta.write(fasta_str + "\n")
            except FileExistsError:
                print(f"{file_path.name} already exists, you can find it in {file_path}")


def is_not_hidden_file(filename:str) -> bool:
    """
    Ignore hidden files when scanning folders.
    """
    return not filename.startswith(".")


def generate_blast_hit_table(fasta_path: Path, blast_results_dir: Path) -> None:
    """
    Run self-BLAST for each FASTA file.
    """
    blast_results_dir.mkdir(parents=True, exist_ok=True)

    hit_table_path = blast_results_dir / f"{fasta_path.stem}.txt"

    if not hit_table_path.exists():
        cmd = f"""
        source /etc/profile
        module load BLAST

        blastn -query "{str(fasta_path)}" \
        -subject "{str(fasta_path)}" \
        -out "{str(hit_table_path)}" \
        -outfmt 6
        """

        subprocess.run(["bash", "-lc", cmd], check=True)

    return hit_table_path
