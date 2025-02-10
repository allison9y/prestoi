import os
import subprocess
import argparse
import requests
import json
import glob
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder
from utils import makedir_if_not_exists, is_dir, is_file
from parsers import deduplicate_stockholm_msa, remove_empty_columns_from_stockholm_msa, convert_stockholm_to_a3m, parse_hhr
from Bio import pairwise2
from itertools import product
from statistics import mode

# Step 1: Parse FASTA File
def parse_fasta(input_fasta):
    """Reads a FASTA file and returns a list of sequences."""
    sequences = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequences.append(record)
    return sequences

# Step 2: Run Jackhmmer for MSA search with default parameters
def run_jackhmmer(input_fasta, output_prefix, uniref90_db):
    """Runs Jackhmmer to search against UniRef90 and generate MSA with default parameters."""
    msa_output = f"{output_prefix}.sto"
    tblout_path = f"{output_prefix}.tblout"

    print("Running jackhmmer.....")

    if os.path.exists(msa_output):
        return msa_output

    cmd = [
        'jackhmmer',
        "-o", "/dev/null",  # Suppress stdout
        "-A", msa_output,   # Save alignment in Stockholm format
        "--noali",          # Do not save the full sequence alignment
        "--F1", "0.0005",   # MSV filter threshold
        "--F2", "0.00005",  # Viterbi filter threshold
        "--F3", "0.0000005", # Forward filter threshold
        "--incE", "0.0001",  # Inclusion E-value
        "-E", "0.0001",      # E-value threshold
        "--cpu", "8",        # Number of CPU cores
        "-N", "1",           # Number of iterations
        "--tblout", tblout_path,  # Save tabular output
        input_fasta,
        uniref90_db
    ]
    
    subprocess.run(cmd, check=True)
    return msa_output

# Step 3: Run HHmake to generate HMM profiles
def run_hhmake(hhmake_binary, msa_file, output_hmm):
    """Runs HHmake to generate an HMM profile from MSA."""
    print("Running hhmake.....")
    if os.path.exists(output_hmm):
        return

    cmd = [hhmake_binary, "-i", msa_file, "-o", output_hmm]
    subprocess.run(cmd, check=True)
    return output_hmm

# Step 4: Run HHsearch to search for templates
def run_hhsearch(hhsearch_binary, hmm_file, output_hhr, hhdb_prefix):
    """Runs HHsearch to find templates using the HMM profile with default parameters."""
    # Ensure the database exists
    if not glob.glob(hhdb_prefix + '_*'):
        raise ValueError(f"Could not find HHsearch database at {hhdb_prefix}")

    print("Running hhsearch.....")

    if not os.path.exists(output_hhr): 
        MAXSEQ = 1_000_000  # Maximum number of sequences in alignment

        cmd = [
            hhsearch_binary,
            "-i", hmm_file,          # Input HMM file
            "-o", output_hhr,        # Output HHR file
            "-maxseq", str(MAXSEQ),  # Maximum sequences in alignment
            "-d", hhdb_prefix        # HHsearch database prefix
        ]

        subprocess.run(cmd, check=True)

    with open(output_hhr) as f:
        hhr = f.read()
    return hhr


def get_pdb_stoichiometry(pdb_codes):
    """
    Retrieves global stoichiometry information for a list of PDB codes from the RCSB PDB API.
    """
    base_url = "https://data.rcsb.org/rest/v1/core/assembly/"
    stoichiometry_info = {}

    for pdb in pdb_codes:
        url = f"{base_url}{pdb}/1"
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()
            stoichiometry = data.get("rcsb_struct_symmetry", [{}])[0].get("stoichiometry", "Not Available")
            if stoichiometry == "Not Available":
                stoichiometry = data.get("pdbx_struct_assembly", {}).get("details", "Not Available")
            stoichiometry_info[pdb] = stoichiometry
        else:
            stoichiometry_info[pdb] = "Not Found"

    return stoichiometry_info

def extract_chain_sequences(pdb_file):
    """
    Extracts sequences of chains from a given PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        dict: Dictionary mapping chain IDs to their sequences.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ppb = PPBuilder()
    chain_sequences = {}

    for model in structure:
        for chain in model:
            peptides = ppb.build_peptides(chain)
            if peptides:  # Check if there are polypeptides in the chain
                sequence = "".join([str(peptide.get_sequence()) for peptide in peptides])
                chain_sequences[chain.id] = sequence

    return chain_sequences

from Bio import pairwise2
from collections import defaultdict

def compute_sequence_identity(seq1, seq2):
    """Compute percentage sequence identity between two sequences using globalxx."""
    # If you prefer scoring-based alignment (globalms), uncomment the below.
    # alignment = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1, one_alignment_only=True)[0]
    # matches = sum(1 for a, b in zip(alignment.seqA, alignment.seqB) if a == b)
    # identity = matches / max(len(seq1), len(seq2))

    alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    identity = alignment.score / max(len(seq1), len(seq2))
    return identity

def infer_stoichiometry_from_chains(chain_sequences, identity_threshold=0.9):
    """
    Infers stoichiometry from chain sequences by clustering similar chains.

    Args:
        chain_sequences (dict): Dictionary mapping chain IDs to sequences.
        identity_threshold (float): Threshold for sequence identity to cluster chains (default 90%).

    Returns:
        (str, dict):
          str: Inferred stoichiometry (e.g., "A5" for five nearly identical chains).
          dict: Mapping from cluster label -> list of chain IDs in that cluster.
    """
    # cluster_representatives maps label -> list_of_chain_ids
    cluster_representatives = {}
    chain_labels = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZ")  # Label generator

    # Convert dict to list to iterate items consistently
    chain_items = list(chain_sequences.items())  # [(chain_id, chain_seq), ...]

    for chain_id, chain_seq in chain_items:
        assigned = False

        # Try to assign chain to an existing cluster
        for label, chain_ids in cluster_representatives.items():
            # representative = chain_ids[0] is the first chain's ID in that cluster
            rep_id = chain_ids[0]
            rep_seq = chain_sequences[rep_id]

            identity = compute_sequence_identity(chain_seq, rep_seq)
            # Debug print if desired:
            # print(f"Comparing: {chain_id} seq vs {rep_id} seq; identity={identity}")

            if identity >= identity_threshold:
                # Add chain_id to the existing cluster
                cluster_representatives[label].append(chain_id)
                assigned = True
                break

        if not assigned:
            # Create a new cluster with the next available label
            new_label = next(chain_labels)
            cluster_representatives[new_label] = [chain_id]

    # Format the stoichiometry output, e.g. A2B3...
    stoichiometry_str = "".join(
        f"{label}{len(chain_ids)}"
        for label, chain_ids in sorted(cluster_representatives.items())
    )

    return stoichiometry_str, cluster_representatives

def match_chains_to_stoichiometry(input_sequence, template_name, pdb_file):
    """
    Groups chain IDs into components based on sequence identity, infers stoichiometry,
    and finds the best matching component for the input sequence.

    Args:
        input_sequence (str): Input sequence to compare.
        pdb_file (str): Path to the PDB file.

    Returns:
        dict: {
          'inferred_stoichiometry': str,
          'component_to_chains': dict,
          'best_component': str,
          'best_identity': float,
          'num_of_copies': int
        }
    """
    # Extract chain sequences from PDB
    chain_sequences = extract_chain_sequences(pdb_file)

    global_stoichiometry = get_pdb_stoichiometry([template_name]).get(template_name)[0]
    if global_stoichiometry and re.match(r"^A\d+$", global_stoichiometry):
        # If stoichiometry indicates a homo-multimer (e.g., A2, A4, etc.)
        parsed_stoichiometry = parse_stoichiometry(global_stoichiometry)
        num_of_copies = parsed_stoichiometry.get("A", 0)  # Get the number of copies for 'A'
        # print(f"{template_name}: Homo-multimer detected with stoichiometry {global_stoichiometry}. "
        #         f"Using {num_of_copies} as the number of copies.")
        return {
            "inferred_stoichiometry": global_stoichiometry,
            "component_to_chains": {"A": list(chain_sequences.keys())},
            "best_component": "A",
            "best_identity": 1.0,  # Assume identity is perfect for homo-multimers
            "num_of_copies": num_of_copies
        }
        
    # Infer stoichiometry and clustering from chain IDs
    inferred_stoichiometry, cluster_map = infer_stoichiometry_from_chains(chain_sequences)

    # cluster_map is like {'A': ['Chain1', 'Chain2'], 'B': ['Chain3', 'Chain4']}
    # We'll rename it here for clarity
    component_to_chains = cluster_map

    # Find best match for input_sequence
    best_component = None
    best_identity = 0.0

    for component_label, chain_ids in component_to_chains.items():
        for chain_id in chain_ids:
            chain_seq = chain_sequences[chain_id]
            identity = compute_sequence_identity(input_sequence, chain_seq)

            if identity > best_identity:
                best_identity = identity
                best_component = component_label

    return {
        "inferred_stoichiometry": inferred_stoichiometry,
        "component_to_chains": component_to_chains,
        "best_component": best_component,
        "best_identity": best_identity,
        "num_of_copies": len(component_to_chains[best_component]) if best_component else 0,
    }

def determine_subunit_copies(subunit_templates, pdb_folder):
    """
    Determines the possible number of copies for each subunit based on sequence-template matches.

    Args:
        subunit_templates (dict): { subunit -> list of template hits }
        pdb_folder (str): Path containing PDB files.

    Returns:
        dict: { subunit -> list of possible copy numbers }
    """
    unique_pdbs = set(template.name[:4].lower() for templates in subunit_templates.values() for template in templates)
    subunit_possible_counts = {}

    os.makedirs(pdb_folder, exist_ok=True)

    for subunit, templates in subunit_templates.items():
        print(f"Listing possible copies of the subunit {subunit}")
        subunit_possible_counts[subunit] = []

        for template in templates:
            template_name = template.name.lower()[:4]
            pdb_file = os.path.join(pdb_folder, f"{template_name}.pdb1")

            # Download missing PDB files if needed
            if not os.path.exists(pdb_file):
                os.system(f"wget -q https://files.rcsb.org/download/{template_name}.pdb1 -P {pdb_folder}")
                if not os.path.exists(pdb_file):
                    # Could not retrieve the file; skip
                    print(f"{template_name}: cannot download the pdb1 for this template")
                    continue

            # Directly match subunit sequence to stoichiometry from the PDB's chain sequences
            match_dict = match_chains_to_stoichiometry(template.hit_sequence.replace('-', ''), template_name, pdb_file)

            copies = match_dict['num_of_copies']
            print(f"{template_name}: stoichiometry: {match_dict['inferred_stoichiometry']}, "
                  f"possible number of copies for {subunit}: {copies}")

            subunit_possible_counts[subunit].append(copies)

    return subunit_possible_counts

def generate_combinations(possible_copies):
    subunit_combinations = [(subunit.split('_')[-1], sorted(set(copies))) for subunit, copies in possible_copies.items()]
    cartesian_product = product(*[[(subunit, count) for count in counts] for subunit, counts in subunit_combinations])
    combinations = [''.join(f"{subunit}{count}" for subunit, count in combination) for combination in cartesian_product]
    return combinations

def parse_stoichiometry(stoichiometry):
    """
    Parses stoichiometry string into a dictionary (e.g., 'A2B3' -> {'A': 2, 'B': 3}).
    """
    stoich_dict = {}
    
    matches = re.findall(r'([A-Z])(\d*)', stoichiometry)

    for subunit, count in matches:
        stoich_dict[subunit] = int(count) if count else 1

    return stoich_dict

def get_single_template_for_all_subunits(subunit_templates):
    """
    Checks if there is one template that appears in *all* subunits' top hits.
    Returns that template name if found, otherwise None.
    """
    all_sets = []
    for templates in subunit_templates.values():
        # Each template in the list is an object with .name
        # Extract the 4-letter code or use the entire .name
        template_ids = {t.name[:4].lower() for t in templates}
        all_sets.append(template_ids)

    # Intersection of all sets
    common = set.intersection(*all_sets)
    return next(iter(common)) if common else None

def partition_copies(total_copies, num_subunits):
    """
    Evenly distributes total copies among the given number of subunits.

    Args:
        total_copies (int): Total number of copies to distribute.
        num_subunits (int): Number of subunits.

    Returns:
        list: A list where each element represents the number of copies assigned to a subunit.
    """
    base_count = total_copies // num_subunits  # Minimum copies each subunit gets
    remainder = total_copies % num_subunits  # Extra copies to distribute

    # Distribute the extra copies among the first 'remainder' subunits
    distribution = [base_count + 1] * remainder + [base_count] * (num_subunits - remainder)

    return distribution

########################################
# SINGLE TEMPLATE COVERAGE ACROSS ALL SUBUNITS
########################################
def get_single_template_for_all_subunits(subunit_templates):
    """
    Finds a single template that appears in *all* subunits' top hits.
    Returns the template with the highest average sum_probs across subunits.

    Args:
        subunit_templates (dict): {subunit_name -> list of template objects}.

    Returns:
        str: Template code with the highest average sum_probs (or None if no common template).
    """
    all_sets = []
    template_scores = defaultdict(list)  # Store sum_probs for each template across subunits

    for templates in subunit_templates.values():
        tset = {t.name[:4].lower() for t in templates}  # Extract 4-letter PDB codes
        all_sets.append(tset)

        # Store sum_probs for each template
        for t in templates:
            template_scores[t.name[:4].lower()].append(t.sum_probs)

    # Find templates that appear in all subunits
    common = set.intersection(*all_sets)
    if not common:
        return None  # No common template found

    # Compute average sum_probs for common templates
    sorted_templates = sorted(
        common, 
        key=lambda x: sum(template_scores[x]) / len(template_scores[x]), 
        reverse=True  # Sort in descending order (highest sum_probs first)
    )
    
    return sorted_templates

def handle_single_template_stoich(
    template_code,
    subunit_list,
    subunit_sequences,
    pdb_folder,
    identity_threshold=0.3
):
    """
    1) Download the PDB if needed.
    2) Extract chain sequences -> infer stoichiometry from cluster_map.
    3) For each cluster label (A, B, etc.), see how many subunits match it.
       - If multiple subunits match the same label, partition the total copies among them.
    4) Build a final distribution string, e.g. "H0208_A:2,H0208_B:2".
    """

    # Ensure PDB file is downloaded
    pdb_file = os.path.join(pdb_folder, f"{template_code}.pdb1")
    if not os.path.exists(pdb_file):
        os.system(f"wget -q https://files.rcsb.org/download/{template_code}.pdb1 -P {pdb_folder}")
        if not os.path.exists(pdb_file):
            print(f"Cannot retrieve {template_code}.pdb1")
            return None

    # 1) Extract chain sequences
    chain_seqs = extract_chain_sequences(pdb_file)
    if not chain_seqs:
        print(f"No chain sequences found for template {template_code}")
        return None

    if len(chain_seqs) < len(subunit_sequences):
        print(f"The chain number in {template_code} is less than the number of input subunits")
        return None

    # 2) Infer stoichiometry from chain IDs -> we get a string (e.g. "A3B2") and cluster_map
    inferred_stoich_str, cluster_map = infer_stoichiometry_from_chains(chain_seqs)
    # Example:
    #   inferred_stoich_str = "A3B2"
    #   cluster_map = {'A': ['Chain1','Chain2','Chain3'], 'B': ['Chain4','Chain5']}

    # Parse that stoich string into a dict, e.g. {'A':3, 'B':2}
    parsed_stoich = parse_stoichiometry(inferred_stoich_str)
    # print(chain_seqs)

    # 3) Identify which cluster label each subunit matches best
    #    We'll compare subunit sequences to each cluster label's chain sequences
    #    and pick the label that yields the highest identity.
    label_to_subunits = defaultdict(list)

    for sub_name in subunit_list:
        sub_seq = subunit_sequences[sub_name]
        best_label = None
        best_identity = 0.0

        for label, chain_ids in cluster_map.items():
            # chain_ids is a list of chain IDs in that cluster
            # We'll find the best identity among them
            local_best = 0.0
            for chain_id in chain_ids:
                cseq = chain_seqs[chain_id]
                ident = compute_sequence_identity(sub_seq, cseq)
                if ident > local_best:
                    local_best = ident
            if local_best > best_identity:
                best_identity = local_best
                best_label = label

        if best_label:
            label_to_subunits[best_label].append(sub_name)

    # print(label_to_subunits)
    # 4) For each cluster label => we have a total # of copies (parsed_stoich[label])
    #    If multiple subunits map to that label, we partition. If 1 subunit => it gets them all.
    final_assignments = {}

    for label, total_copies in parsed_stoich.items():
        # subunits that matched label
        matched_subs = label_to_subunits[label]
        # print(matched_subs)
        if not matched_subs:
            # no subunits matched => skip
            final_assignments[label] = {}
            continue

        if len(matched_subs) == 1:
            # single subunit => it gets all copies
            final_assignments[label] = {matched_subs[0]: total_copies}
        else:
            # partition among multiple subunits
            # e.g. if total_copies=4 and matched_subs=2 => possible combos [1,3],[2,2],[3,1]
            distribution = partition_copies(total_copies, len(matched_subs))

            # Save the best distribution
            final_assignments[label] = {}
            for i, subn in enumerate(matched_subs):
                final_assignments[label][subn] = distribution[i]

    # 5) Build a final distribution string
    #    e.g. "H0208_A:2,H0208_B:2"
    stoich_parts = []
    for label, sub_dict in final_assignments.items():
        for subn, cnt in sub_dict.items():
            stoich_parts.append(f"{subn.split('_')[-1]}{cnt}")
    final_str = "".join(sorted(stoich_parts))

    print(f"[handle_single_template_stoich] Template={template_code}, "
          f"cluster_stoich={inferred_stoich_str}, final_distribution={final_str}")
    return final_str

########################################
# FINALIZE CONFIDENT STOICHIOMETRY
########################################
def finalize_confident_stoichiometry(subunit_templates, possible_copies, pdb_folder, subunit_sequences):
    """
    1) If a single template covers all subunits, parse its stoichiometry.
       - If it's homo or heteromultimer, generate partitions accordingly.
    2) If no single template found, fallback to mode-based approach.
    """
    single_templates = get_single_template_for_all_subunits(subunit_templates)
    if single_templates is not None:
        for single_template in single_templates:
            # Attempt to interpret the single template stoichiometry for all subunits
            # subunit_list = sorted(subunit_templates.keys()) # e.g. ["H0208_A","H0208_B"]
            subunit_list = list(subunit_templates.keys())
            final_str = handle_single_template_stoich(
                template_code=single_template,
                subunit_list=subunit_list,
                subunit_sequences=subunit_sequences,  # We assume you have subunit->seq
                pdb_folder=pdb_folder,
                identity_threshold=0.3
            )
            if final_str is not None and len(final_str) > 0:
                return final_str

    return ""

def process_unique_sequences(sequences, output_path, uniref90_db, hhmake_binary, hhsearch_binary, hhdb_prefix):
    """Process unique sequences for homomultimer optimization."""
    sequence_map = {}
    unique_sequences = {}
    subunit_templates = {}

    for subunit, sequence in sequences.items():
        if sequence in unique_sequences:
            sequence_map[subunit] = unique_sequences[sequence]
        else:
            sequence_map[subunit] = subunit
            unique_sequences[sequence] = subunit

            # Generate files and process sequence
            fasta_path = os.path.join(output_path, f"{subunit}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{subunit}\n{sequence}\n")

            msa_file = run_jackhmmer(fasta_path, os.path.join(output_path, subunit), uniref90_db)
            a3m_file = os.path.join(output_path, f"{subunit}.a3m")
            with open(a3m_file, "w") as f:
                msa_content = open(msa_file).read()
                msa_content = deduplicate_stockholm_msa(msa_content)
                msa_content = remove_empty_columns_from_stockholm_msa(msa_content)
                msa_content = convert_stockholm_to_a3m(msa_content)
                f.write(msa_content)

            hmm_file = os.path.join(output_path, f"{subunit}.hmm")
            run_hhmake(hhmake_binary, a3m_file, hmm_file)

            hhr_file = os.path.join(output_path, f"{subunit}.hhr")
            hhsearch_result = run_hhsearch(hhsearch_binary, hmm_file, hhr_file, hhdb_prefix)
            parsed_templates = parse_hhr(hhsearch_result)

            subunit_templates[subunit] = sorted(parsed_templates, key=lambda x: x.sum_probs, reverse=True)[:10]

    return subunit_templates, sequence_map

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta',type=str,required=True)
    parser.add_argument('--output_path', type=str, required=True)

    args = parser.parse_args()
    ## Get Alphafold3 configurations
    if not os.path.exists("config.json"):
        print("Please run configure_af3.py first to create a config.json file containing af3_program_path,af3_params_path, and af3_db_path")
        exit()
    with open("config.json", "r") as f:
        config = json.load(f)

    os.makedirs(args.output_path, exist_ok=True)
    af3_db_path = config["af3_db_path"]
    uniref90_database = os.path.join(af3_db_path, 'uniref90_2022_05.fa')
    hhdb_prefix = config["hhdb_prefix"]
    hhsearch_binary = config["hhsearch_binary"]
    hhmake_binary = config["hhmake_binary"]

    if not os.path.exists(uniref90_database):
        print("Please run configure.py first to create a config.json file containing af3_program_path,af3_params_path, and af3_db_path")
        exit()

    sequences = parse_fasta(args.input_fasta)

    sequences = parse_fasta(args.input_fasta)
    subunit_sequences = {record.id: str(record.seq) for record in sequences}

    subunit_templates, sequence_map = process_unique_sequences(
        subunit_sequences,
        args.output_path,
        uniref90_database,
        hhmake_binary,
        hhsearch_binary,
        hhdb_prefix
    )
    
    # Determine possible subunit copies based on templates
    pdb_folder = os.path.join(args.output_path, 'templates')
    os.makedirs(pdb_folder, exist_ok=True)
    possible_copies = determine_subunit_copies(subunit_templates, pdb_folder)

    # Print possible copy numbers for each subunit
    # print(possible_copies)
    for subunit, copies in possible_copies.items():
        unique_copies = set(copies)
        print(f"Subunit {subunit}: Possible Copies - {sorted(unique_copies)}")
    
    combos = generate_combinations(possible_copies)
    print("Stoichiometry candidates:", combos)
    
    # Final step: Single-template logic vs. fallback
    final_stoichiometry = finalize_confident_stoichiometry(
        subunit_templates=subunit_templates,
        possible_copies=possible_copies,
        pdb_folder=pdb_folder,
        subunit_sequences=subunit_sequences
    )
    if len(final_stoichiometry) > 0:
        print("Template-based stoichiometry prediction:", final_stoichiometry)
    else:
        print("No sufficient evidence to make template-based stoichiometry prediction")
    


