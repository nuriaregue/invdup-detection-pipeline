"""
This script analyses each read and uses the information stored in SA tags to identify and recover supplementary alignments that are not already captured in the inversion region.

SA entries are first compared with the regional alignments already
captured for the read. The MAPQ filter is then applied only to SA
entries describing alignments that are not already present.

Only new SA entries with MAPQ >= 20 are considered reliable and used
for alignment recovery.

SA recovery statuses:
    - NO_SA:The read does not contain any SA tag.
    - SA_ALREADY_IN_REGION:All SA entries correspond to alignments already captured in the
    inversion region.
    - NO_RELIABLE_SA: The read contains new SA entries, but none passes the MAPQ filter.
    - SA_RECOVERED:At least one new reliable SA alignment was recovered from the BAM.
    - SA_OTHER_CHROMOSOMES:All new reliable SA entries point to chromosomes different from
    the inversion chromosome and were not recovered.
    - SA_SAME_CHR_NOT_RECOVERED:At least one new reliable SA entry points to the inversion
    chromosome, but its alignment was not recovered from the BAM.

Recovered SA alignments are stored separately and do not modify the
original regional pattern classification.
"""


from pathlib import Path
from scripts.bam_analysis import utils_bam as ub
from scripts.config import (
    NCBI_TO_CHR
)
import pysam

def normalize_reference_name(ref_name: str) -> str:
    """
    Normalize reference names to a chr-style format when possible.

    Examples:
        chr1 -> chr1
        1 -> chr1
        X -> chrX
        ref|NC_000023.11|:1-156040895 -> chrX
    """
    if ref_name.startswith("chr"):
        return ref_name

    for ncbi_name, chr_name in NCBI_TO_CHR.items():
        if ncbi_name in ref_name:
            return chr_name

    if ref_name in {"X", "Y"}:
        return f"chr{ref_name}"

    if ref_name.isdigit():
        return f"chr{ref_name}"

    return ref_name

def get_sa_key(sa_dict: dict) -> tuple:
    """
    Build a unique key for one SA entry.
    """
    return (
        sa_dict["chr"],
        sa_dict["start"],
        sa_dict["strand"],
        sa_dict["cigar"],
    )

def get_read_sa_tags(read: dict) -> list:
    """
    Each alignment of the same read can contain an SA tag pointing to the other
    supplementary alignments, so the same SA entry may appear more than once.
    We deduplicate SA entries to avoid counting the same alignment repeatedly.
    """
    sa_tags = {}
    for alignment in read.get("alignments", []):
        if alignment.has_tag("SA"):
            sa_list_dict = sa_tag_to_dict(alignment.get_tag("SA"))
            for sa_dict in sa_list_dict:
                sa_key=get_sa_key(sa_dict)
                sa_tags[sa_key]=sa_dict
    return list(sa_tags.values())

def sa_tag_to_dict(sa_tag:str):
    """
    SA_TAG is formatted as a semicolon-delimited list: (rname, pos, strand, CIGAR, mapQ, NM;)
    It is converted to a dict with the same fields
    """
    sa_list=[]
    aligns = sa_tag.strip(";").split(";")
    for al in aligns:
        fields=al.split(",")
        if len(fields)==6:
            # Convert from ref|NC_000023.11|:1-156040895 to its corresponding chromosome
            # If the format is correct it won't be processed
            normalized_chr = normalize_reference_name(fields[0]) 
            sa_dict={"chr":normalized_chr,"start":int(fields[1]),"strand":fields[2], "cigar":fields[3],"mapq":int(fields[4]), "NM":int(fields[5])}
            sa_list.append(sa_dict)
        else:
            print(f"Malformed SA: {al}")
        
    return sa_list

def sa_is_already_in_region(sa: dict, read: dict) -> bool:
    found = False

    for alignment in read["alignments"]:
        if alignment.reference_name == sa["chr"] and alignment.reference_start == sa["start"] - 1:
            found = True

    return found


def search_supplementary(classified_reads:dict,bam_path:Path,inv_chr:str):
    """
    Method that searches for supplementary alignments and creates read["recovered_sa_alignments"] fields
    """
    
    for read_id, read in classified_reads.items():
        read["unrecovered_sa_tags"] = []
        read["n_recovered_sa_alignments"]=0
        read["recovered_sa_alignments"] = []
        read["recovered_sa_flags"] = []
        #Check whether the original regional alignments contain any SA tag,independently of its MAPQ.
        has_sa_tag = False
        
        for alignment in read.get("alignments", []):
            if alignment.has_tag("SA"):
                has_sa_tag = True

        #Get all deduplicated SA tags 
        sa_tags = get_read_sa_tags(read)
        #Keep only SA tags describing alignments not already captured in the region.
        new_sa_tags = []

        for sa in sa_tags:
            if not sa_is_already_in_region(sa, read):
                new_sa_tags.append(sa)
        # From the new SA tags, keep only reliable entries.
        sa_to_recover = []

        for sa in new_sa_tags:
            if sa["mapq"] >= 20:
                sa_to_recover.append(sa)
                #If the record is not present in the subset BAM, we keep the SA tag information so it can later be added to segments.tsv as SA_TAG_ONLY evidence.
                record_found = fetch_supplementaries(sa, bam_path, read_id, classified_reads)
                if not record_found:
                    read["unrecovered_sa_tags"].append(sa)
        read["sa_tags"] = sa_tags
        read["n_sa_tags"] = len(sa_tags)
        read["has_sa_tag"] = has_sa_tag
        read["has_reliable_sa"] = len(sa_to_recover) > 0
        read["new_sa_tags"] = new_sa_tags
        read["n_new_sa_tags"] = len(new_sa_tags)
        read["sa_chromosomes"] = sorted({sa["chr"] for sa in sa_to_recover})
        read["n_sa_to_recover"] = len(sa_to_recover)
        read["n_unrecovered_sa_tags"] = len(read["unrecovered_sa_tags"])
       
        #Only reliable SA tags pointing outside the inversion region are considered (mapq>=20)
        if not read["has_sa_tag"]:  #The read does not contain any SA tag.
            read["sa_recovery_status"] = "NO_SA"
        elif read["n_new_sa_tags"] == 0: #All SA entries correspond to alignments already captured in the region.
            read["sa_recovery_status"] = "SA_ALREADY_IN_REGION"
        elif not read["has_reliable_sa"]: #The read has SA_tags, but none passes the MAPQ filter we have defined.
            read["sa_recovery_status"] = "NO_RELIABLE_SA"
        elif read["n_recovered_sa_alignments"] > 0: #New alignments were found
            read["sa_recovery_status"] = "SA_RECOVERED"
        else:
            other_chrm = only_other_chromosomes(sa_to_recover, inv_chr)
            if other_chrm: 
                #The reliable SA alignments point only to other chromosomes.
                read["sa_recovery_status"] = "SA_OTHER_CHROMOSOMES"
            else:
                #A reliable SA alignment on the inversion chromosome was not recovered.
                read["sa_recovery_status"] = "SA_SAME_CHR_NOT_RECOVERED"

def get_alignment_key(alignment) -> tuple:
    """
    Build a comparable key for one alignment already captured from the BAM.
    """
    strand = "-" if alignment.is_reverse else "+"

    return (
        alignment.reference_name,
        alignment.reference_start,
        strand,
        alignment.cigarstring,
    )

def fetch_supplementaries(sa_tag:dict,bam_path:Path,read_id:str,classified_reads:dict):
    margin=1000
    record_found = False
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        chromosome=sa_tag["chr"]
        region_start= sa_tag["start"] -1 #start is 1-based so we convert it to 0-based
        window_start= max(0,region_start-margin)
        window_end=region_start+margin
        if chromosome in bamfile.references:
            alignments_iter = bamfile.fetch(chromosome, window_start, window_end)
            for align in alignments_iter:
                #To confirm that the alignment found in the fetch window matches the SA tag
                # we compare its read_id, strand and its position .
                query_name = align.query_name
                align_strand = "-" if align.is_reverse else "+"
                same_position = abs(align.reference_start - region_start) <= 5
                same_strand = align_strand == sa_tag["strand"]
                if query_name == read_id and same_position and same_strand:
                    add_recovered_sa_alignment(classified_reads[read_id],align)
                    record_found = True

        else:
            print(f"Reference not found in BAM: {chromosome}")

    return record_found

def add_recovered_sa_alignment(read: dict, alignment) -> bool:
    """
    Save an SA-recovered alignment separately from the original regional
    alignments. This keeps the original flag-based classification unchanged.
    """

    existing_keys = []

    for existing_align in read["alignments"]:
        existing_keys.append(get_alignment_key(existing_align))

    for existing_align in read["recovered_sa_alignments"]:
        existing_keys.append(get_alignment_key(existing_align))

    alignment_key = get_alignment_key(alignment)

    added = False

    if alignment_key not in existing_keys:
        read["recovered_sa_alignments"].append(alignment)
        read["recovered_sa_flags"].append(alignment.flag)
        read["n_recovered_sa_alignments"] = len(read["recovered_sa_alignments"])
        added = True

    return added


def print_read(classified_read: dict) -> str:
    alignments=classified_read.get("alignments",[])
    if len(alignments) == 0:
        print("\tRead without alignments")
    else: 
        read_id = alignments[0].query_name
        print(f"\t {read_id}|{classified_read['n_alignments']} alignments |" 
            f"flags = {classified_read['flags']}|pattern = {classified_read['pattern']}|")
    
def only_other_chromosomes(sa_tags: list, inv_chr: str) -> bool:
    """
    Return True if all SA tags point to chromosomes different from the inversion chromosome.
    """
    for sa in sa_tags:
        if sa["chr"] == inv_chr:
            return False
    return True
