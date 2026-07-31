"""
Central configuration for the INV_DUP analysis project.

It contains shared project paths and parameters used by the BAM
analysis, read geometry analysis and self-BLAST pipelines.
"""

from pathlib import Path


# scripts/config.py -> project root is one folder above scripts/
PROJECT_DIR = Path(__file__).resolve().parents[1]


# Main directories
BAM_DIR = PROJECT_DIR / "BAM"
REF_DIR = PROJECT_DIR / "REF"

PIPELINE_INV_DUP_RESULTS_DIR = PROJECT_DIR / "results" / "invdup_filtering"
PIPELINE_SELF_BLAST_RESULTS_DIR = PROJECT_DIR / "results" / "self_blast"

#Self blast directories
FASTA_DIR = PIPELINE_SELF_BLAST_RESULTS_DIR / "FASTA"
BLAST_RESULTS_DIR = PIPELINE_SELF_BLAST_RESULTS_DIR / "BLAST_RESULTS"
OUTPUT_TSV_DIR = PIPELINE_SELF_BLAST_RESULTS_DIR / "OUTPUT_TSV"
BED_DIR = PIPELINE_SELF_BLAST_RESULTS_DIR / "BED"

# Input files
ALLCOORDS_PATH = BAM_DIR / "allcoords.txt"
REFERENCE_PATH = REF_DIR / "Std_ref2.fa"


# BAM analysis outputs used when scripts are run individually
SUMMARY_TSV = BAM_DIR / "read_summary.tsv"
SEGMENTS_TSV = BAM_DIR / "read_segments.tsv"
GEOM_TSV = BAM_DIR / "segments_geometry.tsv"


# Self-BLAST outputs
BLAST_HITS_TSV = OUTPUT_TSV_DIR / "blast_hits.tsv"
BLAST_PAIRS_TSV = OUTPUT_TSV_DIR / "blast_pairs.tsv"
BLAST_SUMMARY_TSV = OUTPUT_TSV_DIR / "blast_summary.tsv"

# Shared parameters
MARGIN = 50000

# Self-BLAST parameters
MIN_IDENTITY = 90.0
MIN_LENGTH = 1000


# Reference-name conversion
NCBI_TO_CHR = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
}
MIN_GEOMETRY_MAPQ = 55
MIN_SA_RECOVERY_MAPQ = 20
MIN_GAP_COVERAGE_FRACTION = 0.80
COORDINATE_MARGIN = 50
MIN_READS_INV_DUP=3