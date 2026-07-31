import pysam
from pathlib import Path

from scripts.config import (
    MIN_GEOMETRY_MAPQ
)


def proces_allcoords(allcoords_path:Path,margin:int)->dict:
    all_inv = {}
    # For each line in allcoords, we need to find the inversion line
    with allcoords_path.open("r") as file_inv:
        #We don't need to store the first line
        header = file_inv.readline().strip().split("\t")
        for line in file_inv:
            #If the line is not empty, we extract the inversion information
            if line.strip(): 
                columns = line.strip().split("\t")
                #Once we have each line of information, we get our fields of interest
                inv_id = columns[0]
                chromosome = "chr" + columns[1]
                #allcoords uses 1-based coordinates, while pysam fetch uses a 0-based start.
                start = int(columns[2]) - 1 
                start=start-margin
                end = int(columns[5])
                end=end + margin
                # For each inversion, we store the chromosome, the start and the end coordinate
                all_inv[inv_id] = [chromosome, start, end]
    return all_inv

def process_reads(bam_path:Path, allcoords_inv:dict,margin:int)->dict:
    #We check that the BAM and inversion coordinates are available.
    if not bam_path.exists() :
        raise FileNotFoundError(f" {bam_path} was not found")
    elif len(allcoords_inv)==0:
        raise KeyError(f"{allcoords_inv} is empty. Execute process_allcoords(ALLCOORDS_PATH) before")

    # Dictionary to store the reads found for each possible inversion
    reads_inv = {}
    classified_reads={}
    # Check if BAM index exists. If not, create it before opening the BAM.
    check_index(bam_path)

    # Open BAM file
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        # We obtain the inversion information that has been stored in the allcoords dictionary 
        for inversion in allcoords_inv:
            inv_info = allcoords_inv.get(inversion)
            chromosome = inv_info[0]
            region_start = int(inv_info[1])
            region_end = int(inv_info[2])
            #We use a dictionary to group alignments by read name
            reads_inv[inversion] = {}
            classified_reads[inversion]={}

            #We get all the alignments that overlap with the inversion region
            try:

                alignments_iter = bamfile.fetch(chromosome, region_start, region_end)

                for read in alignments_iter:
                    query_name = read.query_name

                    #If this read has not been seen before for this possible inversion, 
                    # we create an empty list for its alignments
                    if query_name not in reads_inv[inversion]:
                        reads_inv[inversion][query_name] = []

                    #Store the alignment under its read name
                    reads_inv[inversion][query_name].append(read)

                #Print all reads of a possible inversion
                process_read_align(reads_inv=reads_inv,classified_reads=classified_reads,inversion=inversion)
                    
            except ValueError as err:
                print("Searching for the allcoords.txt inversion line...")
    return classified_reads

def check_index(bam_path:Path):
    """
    Check if BAM index exists. If not, create it before opening the BAM.
    """
    index_name=f"{bam_path}.bai"
    index_bam= bam_path.resolve().parents[0] / index_name
    if not index_bam.exists():
        print(f"BAM index not found for {bam_path}! Creating index...")
        pysam.index(str(bam_path))

def process_read_align(reads_inv:dict,classified_reads:dict,inversion:str):
    if len(reads_inv[inversion]):
        for query_name in reads_inv[inversion]:
            alignments = reads_inv[inversion][query_name]
            alignment_flags = []
                        
            for alignment in alignments:
            #alignment is an AlignedSegment object, so we can access its SAM/BAM fields directly
                alignment_flags.append(alignment.flag) 
                        
                if query_name not in classified_reads[inversion]:
                    reads_inv[inversion][query_name] = {}

                classified_reads[inversion][query_name] = {
                    "n_alignments": len(alignment_flags),
                    "flags": alignment_flags,
                    "pattern": classify_flags(alignment_flags),
                    "alignments": alignments
                }
def classify_flags(flags:list)-> str:
    """
    Given the list of SAM flags for all alignments of the same read,
    return a preliminary pattern classification
    """
    non_dupl_flags=set(flags)
    if len(non_dupl_flags)==1: #[0] or [16] (forward or reverse)
        if 0 in non_dupl_flags or 16 in non_dupl_flags:
            return "SIMPLE"
        elif 2064 in non_dupl_flags or 2048 in non_dupl_flags: #they will be processed later
            return "ONLY_SUPPL"
        else:
            return "UNKNOWN"

    elif len(non_dupl_flags)==2:
        #there is a split, but both alignments have the same direction
        if (16 in non_dupl_flags and 2064 in non_dupl_flags) or (0 in non_dupl_flags and 2048 in non_dupl_flags):
            return "SPLIT"
        #there is a split, but the alignments have different directions
        elif (0 in non_dupl_flags and 2064 in non_dupl_flags) or (16 in non_dupl_flags and 2048 in non_dupl_flags):
            return "SPLIT_INV"
        else:
            return "UNKNOWN"
    else: #Only if we have [0,2048,2064] (len>3) it is considered directy a POSS_INV_DUP
        return "POSS_INV_DUP"
    
def get_reads_by_pattern(reads_by_inv:dict,pattern:str)->dict:
    interesting_reads= {}
    for read_id,read in reads_by_inv.items():
        if read["pattern"] == pattern:
            interesting_reads[read_id]=read
    return interesting_reads

def get_inv_from_bam_name(bam_path: Path) -> str:
    """
    Extract inversion ID from BAM filename

    Example:
    HsInv0036_Std_Inv.bam -> HsInv0036
    HsInv1865_HG00290.bam -> HsInv1865
    """
    return bam_path.stem.split("_")[0]

def classify_alignments_by_mapq_min(alignments:list):
    low_mapq=[]
    min_mapq=[]
    for align in alignments:
        if align.mapq >= MIN_GEOMETRY_MAPQ:
            min_mapq.append(align)
        else:
            low_mapq.append(align)
    return low_mapq,min_mapq

def parse_cigar(cigar: str) -> list:
    """
    Convert a CIGAR string into a list of (length, operation) tuples.

    Example:
        100H500M20I -> [(100, "H"), (500, "M"), (20, "I")]
    """
    operations = []
    number = ""

    for char in cigar:
        if char.isdigit():
            number += char
        else:
            operations.append((int(number), char))
            number = ""

    return operations


def get_reference_length_from_cigar(cigar: str) -> int:
    """
    Return the number of reference bases consumed by the alignment.
    """
    reference_length = 0
    operations = parse_cigar(cigar)

    for length, operation in operations:
        if operation in {"M", "D", "N", "=", "X"}:
            reference_length += length

    return reference_length


def get_query_alignment_length_from_cigar(cigar: str) -> int:
    """
    Return the number of aligned query bases.
    """
    query_alignment_length = 0
    operations = parse_cigar(cigar)

    for length, operation in operations:
        if operation in {"M", "I", "=", "X"}:
            query_alignment_length += length

    return query_alignment_length


def get_query_length_from_cigar(cigar: str) -> int:
    """
    Estimate the complete read length from the CIGAR.
    """
    query_length = 0
    operations = parse_cigar(cigar)

    for length, operation in operations:
        if operation in {"M", "I", "S", "H", "=", "X"}:
            query_length += length

    return query_length

def get_cigar_clipping(cigar: str) -> tuple:
    """
    Return the clipping found at the beginning and end of the CIGAR.
    """
    operations = parse_cigar(cigar)
    leading_clip = 0
    trailing_clip = 0
    start_index = 0

    while start_index < len(operations) and operations[start_index][1] in {"S", "H"}:
        leading_clip += operations[start_index][0]
        start_index += 1
    end_index = len(operations) - 1

    while end_index >= 0 and operations[end_index][1] in {"S", "H"}:
        trailing_clip += operations[end_index][0]
        end_index -= 1
    return leading_clip, trailing_clip
def get_sa_query_coordinates(cigar: str, strand: str) -> dict:
    """
    Infer the SA alignment coordinates in the original read orientation.
    """
    leading_clip, trailing_clip = get_cigar_clipping(cigar)
    alignment_length = get_query_alignment_length_from_cigar(cigar)

    if strand == "+":
        query_start = leading_clip
    elif strand == "-":
        query_start = trailing_clip
    else:
        raise ValueError(f"Invalid strand '{strand}'. Expected '+' or '-'.")

    query_end = query_start + alignment_length

    return {
        "start": query_start,
        "end": query_end
    }