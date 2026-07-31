import csv
from pathlib import Path
from scripts.config import MARGIN
#results_by_invs té una llista de diccionaris hits (hits es la clau, dins de la inv)  
#exemple: results_by_invs[inv001][hits] i aqui cal fer un for per cada hit en hits


def get_hit_intervals(hit: dict):
    """
    Return normalized query and subject intervals from one BLAST hit.
    """
    s_interval={}
    sstart=hit.get("sstart")
    send=hit.get("send")
    if sstart is not None and send is not None:
        s_interval["start"]=min(sstart,send)
        s_interval["end"]=max(sstart,send)
    q_interval={}
    qstart=hit.get("qstart")
    qend=hit.get("qend")
    if qstart is not None and qend is not None:
        q_interval["start"]=min(qstart,qend)
        q_interval["end"]=max(qstart,qend)
    return s_interval,q_interval



def get_hit_pair_key(hit: dict) -> tuple:
    """
    Create a normalized key for one BLAST hit so that reciprocal self-BLAST hits
    are treated as the same repeated block pair.

        1. Extract the query and subject intervals and normalize each interval so start <= end with get_hit_intervals(hit).
        3. Sort the two intervals so A-B and B-A have the same order.
        4. Add the hit orientation to avoid mixing direct and inverted pairs.
    """
    
    s_interval,q_interval = get_hit_intervals(hit) 
    q=(q_interval.get("start"),q_interval.get("end"))
    s = (s_interval.get("start"),s_interval.get("end"))
    sorted_coord=tuple((q,s))
    sorted_coord=tuple(sorted(sorted_coord))
    hit_identification = [sorted_coord,hit["orientation"]]
    sorted_pos = tuple(hit_identification)
    return sorted_pos


def group_blast_hits(hits: list):
    """
    Group raw filtered BLAST hits into unique repeat pairs using exact interval keys.

    Returns a dictionary with for a hit list of a inversion:
    {((qstart1,qend1),(sstart1,send1),orientation):[hit1,hit2],
    ((qstart2,qend2),(sstart2,send2),orientation):[hit3,hit4]}
    """
    unique_hits= {}
    for hit in hits:
        sorted_pos=get_hit_pair_key(hit)
        if is_unique_hit(sorted_pos=sorted_pos,unique_hits=unique_hits):
            unique_hits[sorted_pos]=[]
            unique_hits[sorted_pos].append(hit)
        else:
            unique_hits[sorted_pos].append(hit)
    return unique_hits

def is_unique_hit(sorted_pos:tuple, unique_hits:dict):
    """
    A new unique_hit has a new dictionary entrance.
    TODO: change this method to accept similar hits
    """
    return unique_hits.get(sorted_pos) is None


def convert_to_genomic_coordinates(inv_data:dict, allcoords_inv:dict,inv:str):
    info_allcoords=allcoords_inv.get(inv)
    if info_allcoords is not None:
        pairs=inv_data["pairs"]
        # BLAST positions are local to the extracted FASTA region.
        # Offset (blast start coordinate - 1) converts them back to genomic coordinates.
        offset=info_allcoords["start"] - MARGIN - 1
        if offset<0:
            offset=0
        chrom=info_allcoords["chr"]
        for pair in pairs.values():
            pair["chr"] = chrom
            block_a = pair["block_a"]
            block_b = pair["block_b"]

            pair["block_a"] = (
                block_a[0] + offset,
                block_a[1] + offset,
            )

            pair["block_b"] = (
                block_b[0] + offset,
                block_b[1] + offset,
            )
            rep_hit=pair.get("representative_hit")
            if rep_hit is not None:
                rep_hit["qstart"]=rep_hit["qstart"]+offset
                rep_hit["sstart"]=rep_hit["sstart"]+offset
                rep_hit["qend"]=rep_hit["qend"]+offset
                rep_hit["send"]=rep_hit["send"]+offset
                rep_hit["chr"]=chrom

    
def summarize_blast_hits(hits: list) -> dict:
    """
    Summarize filtered BLAST hits.
    """
    n_direct = 0
    n_inverted = 0
    max_hit_length = 0
    max_hit_identity = 0.0

    for hit in hits:
        if hit["orientation"] == "direct":
            n_direct += 1
        elif hit["orientation"] == "inverted":
            n_inverted += 1

        if hit["length"] > max_hit_length:
            max_hit_length = hit["length"]

        if hit["pident"] > max_hit_identity:
            max_hit_identity = hit["pident"]
    
    summary = {
            "n_hits_filtered": len(hits),
            "n_direct_hits": n_direct,
            "n_inverted_hits": n_inverted,
            "max_hit_length": max_hit_length,
            "max_hit_identity": max_hit_identity,
        }

    return summary

def summarize_blast_paired(pairs:dict):
    """
    paired_hits is [{PAIR_1},{PAIR_2}]
    """

    n_direct = 0
    n_inverted = 0
    max_hit_length = 0
    max_hit_identity = 0.0
    for pair in pairs.values(): 
        repre_hit=pair.get("representative_hit")
        if repre_hit is not None:
            if pair["orientation"] == "direct":
                n_direct += 1
            elif pair["orientation"] == "inverted":
                n_inverted += 1

            if repre_hit["length"] > max_hit_length:
                max_hit_length = repre_hit["length"]

            if repre_hit["pident"] > max_hit_identity:
                max_hit_identity = repre_hit["pident"]
    
    summary = {
            "n_unique_pairs": n_direct+n_inverted,
            "n_unique_direct_pairs": n_direct,
            "n_unique_inverted_pairs": n_inverted,
            "max_pair_hit_length": max_hit_length,
            "max_pair_hit_identity": max_hit_identity,
        }

    return summary

def representative_to_row(pair:dict):
    """
    Convert one deduplicated pair into one TSV row.
    """
    rep_hit=pair.get("representative_hit")
    block_a=pair["block_a"]
    block_b=pair["block_b"]
    if rep_hit is not None:
        inv=pair["inv"]
        row = {
            "PAIR_ID": pair.get("pair_id"),
            "INVERSION":inv,
            "PIDENT": rep_hit["pident"],
            "LENGTH": rep_hit["length"],
            "BLOCK_A_START": block_a[0],
            "BLOCK_A_END": block_a[1],
            "BLOCK_B_START": block_b[0],
            "BLOCK_B_END": block_b[1],
            "ORIENTATION": pair["orientation"],
            "CHROMOSOME":pair["chr"]
        }
        return row

def get_blocks_from_pair_key(pair_key: tuple):
    """
    Extract block coordinates and orientation from the current exact deduplication key.
    """

    blocks, orientation = pair_key
    block_a, block_b = blocks

    return block_a, block_b, orientation

def generate_pair_id(results_by_inv,inv):
    counter=results_by_inv[inv]["pair_counter"]
    pair_id=f"{inv.upper()}_PAIR_{counter}"
    results_by_inv[inv]["pair_counter"]=results_by_inv[inv]["pair_counter"]+1
    return pair_id


def create_pairs_from_groups(grouped_hits:dict,inv:str,results_by_inv:dict):
    """
    Input: a dict with hits grouped by coordinates
    Obtain from grouped hits a pair structure
    """
    """
    {(((1439, 3754), (8544, 10884)), 'inverted'): [{'pident': 95.816, 'length': 2342, 'qstart': 8544, 'qend': 10884, 'sstart': 3754, 'send': 1439, 'orientation': 'inverted'}, 
    {'pident': 95.816, 'length': 2342, 'qstart': 1439, 'qend': 3754, 'sstart': 10884, 'send': 8544, 'orientation': 'inverted'}]}
    """
    pairs = {}
    for pair_key, group in grouped_hits.items():
        block_a, block_b, orientation = get_blocks_from_pair_key(pair_key)

        pair_id = generate_pair_id(results_by_inv, inv)

        group.sort(key=lambda hit: (hit["pident"], hit["length"]), reverse=True)
        representative_hit = group[0]

        pairs[pair_id] = {
            "inv": inv,
            "pair_id": pair_id,
            "block_a": block_a,
            "block_b": block_b,
            "orientation": orientation,
            "raw_hits": group,
            "representative_hit": representative_hit,
        }

    return pairs

def write_blast_pairs_tsv(results_by_inv: dict, output_path: Path,allcoords_dict:dict):
    """
    Save deduplicated repeat pairs for all inversions into a TSV file.
    """
    fieldnames = [
        "PAIR_ID",
        "INVERSION",
        "PIDENT",
        "LENGTH",
        "BLOCK_A_START",
        "BLOCK_A_END",
        "BLOCK_B_START",
        "BLOCK_B_END",
        "ORIENTATION",
        "CHROMOSOME"
        ]
    
    rows = []
    
    for inv in sorted(results_by_inv):

        hits = results_by_inv[inv].get("hits",[])
        unique_hits_structure=group_blast_hits(hits)
        pairs=create_pairs_from_groups(unique_hits_structure,inv,results_by_inv)
        #we add to our global dictionary all the information generated by the deduplication
        results_by_inv[inv]["pairs"]=pairs
        results_by_inv[inv]["summary_pairs"]=summarize_blast_paired(pairs)
        convert_to_genomic_coordinates(results_by_inv[inv],allcoords_dict,inv)
        for pair in pairs.values():
            row = representative_to_row(pair)
            if row is not None:
                rows.append(row)

    with output_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
            delimiter="\t")

        writer.writeheader()
        writer.writerows(rows)
