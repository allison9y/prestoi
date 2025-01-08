import os
import re
import json
from itertools import product

alphabets = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
alphabets = list(alphabets) + [''.join(pair) for pair in product(alphabets, repeat=2)]

def read_fasta(file_path):
    target_name = os.path.splitext(os.path.basename(file_path))[0].lower()
    sequences = []
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                continue
            sequence = line.strip()
            if sequence not in sequences:
                sequences.append(sequence)
            else:
                print("Contains duplicates")
    return target_name, sequences


def create_default_stoichiometry(chain_ids):
    if len(chain_ids) == 1:
        return f"{chain_ids[0]}2"
    default_stoichiometry = "".join([f"{chain}1" for chain in chain_ids])
    return default_stoichiometry

def process_stoichiometry(stoichiometries):
    """Validate if the provided stoichiometries are consistend and extract unique chain IDs"""
    chain_list = []
    for stoichiometry in stoichiometries:
        if not bool(re.fullmatch(r'([a-z]\d+)+', stoichiometry)):
            print("Invalid Stoichiometry Format. Make sure to provide A2,A3... or A1B1,A2B5..etc. (Alphabet followed by a number for each chain)")
            return []
        chains = ','.join(re.findall(r'[A-Za-z]', stoichiometry)).split(',')
        chain_list.append(chains)
    if all(x == chain_list[0] for x in chain_list):
        return chain_list[0]
      
    else:
        print("Inconsistenies in the stoichiometries provided.\nMake sure to provide the same chain IDs with different number of chains.\nFor example:\nAcceptable : A2,A3 or A1B1,A9B18\nUnacceptable: A2,A2B3 or A1B1,A2B2C2")
        return []
    
def validate_sequence_and_stoichiometry(sequences,unique_chains):
    """Validate if there are enough sequences for given stoichiometries """
    if len(sequences)==len(unique_chains):
        return(True)
    else:
        print("Mismatch in number of sequences provided and the given stoichiometries (check both)")
        return(False)

def read_common_data(common_data_path):
    common_data_dict = {}
    with open(common_data_path,"r") as f:
        common_data = json.load(f)
    for chain_info in common_data["sequences"]:
        common_data_dict[chain_info["protein"]["sequence"]] = {}
        common_data_dict[chain_info["protein"]["sequence"]]["unpairedMsa"] = chain_info["protein"]["unpairedMsa"]
        common_data_dict[chain_info["protein"]["sequence"]]["pairedMsa"] = chain_info["protein"]["pairedMsa"]
        common_data_dict[chain_info["protein"]["sequence"]]["templates"] = chain_info["protein"]["templates"]
    
    return common_data_dict 


def generate_json( target_name, stoichiometry, num_seeds, ip_json_path,sequences,common_data,is_default_stoichiometry):

    used_alphabets = 0 
    json_skeleton = {
        "name": f"{target_name}_{stoichiometry}",
        "sequences": [],
        "modelSeeds": [i for i in range(num_seeds)],
        "dialect": "alphafold3",
        "version": 1
    }
    groups = [(group[0], int(group[1])) for group in re.findall(r"([a-z])(\d+)", stoichiometry)]

    # Generate the sequences for the given stoichiometry
    for group_index, (group_letter, count) in enumerate(groups):
        if len(sequences) <= group_index:
            raise ValueError(f"Not enough sequences provided for group {group_letter} in stoichiometry {stoichiometry}.")
            
        sequence = sequences[group_index]  # Get the sequence for the current group
        ids = [alphabets[used_alphabets + i] for i in range(count)]
        
        if is_default_stoichiometry:
            json_skeleton["sequences"].append(
                {
                    "protein": {
                        "id": ids,
                        "sequence": sequence
                    }
                }
            )
        else:
            json_skeleton["sequences"].append(
                {
                    "protein": {
                        "id": ids,
                        "sequence": sequence,
                        "unpairedMsa":common_data[sequence]["unpairedMsa"],
                        "pairedMsa":common_data[sequence]["pairedMsa"],
                        "templates":common_data[sequence]["templates"],
                        
                    }
                }
            )
            

        used_alphabets += count
    json_path = f'{ip_json_path}/{target_name}_{stoichiometry}.json'
    with open(json_path, 'w') as f:
        json.dump(json_skeleton, f,indent=4)
    print("Created : ",json_path)

    return json_skeleton