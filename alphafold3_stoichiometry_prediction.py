import copy
import json
import subprocess
import re
import os
import time
import pandas as pd
import csv
import argparse
from itertools import product
from utils import makedir_if_not_exists,is_dir,is_file
from protein_utils import read_fasta,create_default_stoichiometry,process_stoichiometry,validate_sequence_and_stoichiometry,read_common_data,generate_json


def run_docker(af3_program_path, af3_params_path, af3_db_path, target_name, target_output_path):
    original_cwd = os.getcwd()
    try:
        os.chdir(af3_program_path)
        docker_command = [
            "docker", "run", "--rm",
            "--volume", f"{target_output_path}:/root/af_output",
            "--volume", f"{af3_params_path}:/root/models",
            "--volume", f"{af3_db_path}:/root/public_databases",
            "--gpus", "all",
            "alphafold3",
            "python", "run_alphafold.py",
            f"--json_path=/root/af_output/input_jsons/{target_name}.json",
            "--model_dir=/root/models",
            "--output_dir=/root/af_output"
        ]
        print("Running Docker Command: ", " ".join(docker_command))
        
        try:
            process = subprocess.Popen(
                docker_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1 
            )
            
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print(output.strip())
            
            return_code = process.poll()
            if return_code != 0:
                print(f"Docker command failed with return code {return_code}")
            else:
                print("Docker command completed successfully.")
                
        except Exception as e:
            print(f"Error occurred: {str(e)}")
            
    finally:
        os.chdir(original_cwd)

    
def generate_structures(af3_program_path,af3_params_path,af3_db_path,target_name,target_output_path,sequences,stoichiometries,default_stoichiometry,num_seeds):
    makedir_if_not_exists(target_output_path)
    ip_json_path = os.path.join(target_output_path,"input_jsons")
    makedir_if_not_exists(ip_json_path)
    print(f"####################### Creating input.json for default stoichiometry #######################\n")
    generate_json(target_name,default_stoichiometry,num_seeds,ip_json_path,sequences,common_data={},is_default_stoichiometry=True)

    print(f"####################### Running for default Stoichiometry: {default_stoichiometry} #######################\n")
    print(f"{target_output_path}/{target_name}_{default_stoichiometry}/ranking_scores.csv")
    if not os.path.exists(f"{target_output_path}/{target_name}_{default_stoichiometry}/ranking_scores.csv"):
        run_docker(af3_program_path,af3_params_path,af3_db_path,f"{target_name}_{default_stoichiometry}",target_output_path)
    
    print(f"####################### Completed for Stoichiometry: {default_stoichiometry} #######################\n")

    common_data_path = f"{target_output_path}/{target_name}_{default_stoichiometry}/{target_name}_{default_stoichiometry}_data.json"
    if not os.path.exists(common_data_path):
        print("program Completed without generating results")
        print(common_data_path)
        exit()

    # Processing for other stoichiometries A3,A4......
    print(f"####################### Creating input.json for other stoichiometries #######################\n")

    common_data_dict = read_common_data(common_data_path)
    
    for stoichiometry in stoichiometries:
        if stoichiometry != default_stoichiometry:
            generate_json(target_name,stoichiometry,num_seeds,ip_json_path,sequences,common_data=common_data_dict,is_default_stoichiometry=False)
        
    for stoichiometry in stoichiometries:
        if stoichiometry != default_stoichiometry:
            print(f"####################### Running for Stoichiometry: {stoichiometry} #######################\n")
            ranking_csv_path = f"{target_output_path}/{target_name}_{str(stoichiometry)}/ranking_scores.csv"
            if not os.path.exists(ranking_csv_path):
                run_docker(af3_program_path,af3_params_path,af3_db_path,f"{target_name}_{stoichiometry}",target_output_path)
            print(f"####################### Completed for Stoichiometry: {stoichiometry} #######################\n")


def print_result(target_name,target_output_path,stoichiometries):
    print("\n\nStoichiometry results for : ",target_name.upper())
    print("\nStoichiometry, Maximum ranking score, Average ranking score, Number of models") 
    data = [["stoic","max_ranking_score","avg_ranking_score","num_models"]]   
    for stoichiometry in stoichiometries:
        ranking_csv_path = f"{target_output_path}/{target_name}_{str(stoichiometry)}/ranking_scores.csv"
        if os.path.exists(ranking_csv_path):
            df = pd.read_csv(ranking_csv_path)
            ranking_score_list = df["ranking_score"].tolist()
            max_ranking_score = max(ranking_score_list)
            avg_ranking_score = sum(ranking_score_list)/len(ranking_score_list)
            model_count = len(ranking_score_list)
            print(f"{stoichiometry.upper()},{max_ranking_score},{avg_ranking_score},{model_count}")
            data.append([stoichiometry.upper(),max_ranking_score,avg_ranking_score,model_count])

    result_csv_path = os.path.join(target_output_path,'stoichiometry_results.csv')
    with open(result_csv_path,'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)

    df = pd.read_csv(result_csv_path)

    max_ranking_row = df.loc[df['max_ranking_score'].idxmax()]
    avg_ranking_row = df.loc[df['avg_ranking_score'].idxmax()]
    stoic_max_ranking = max_ranking_row['stoic']
    stoic_avg_ranking = avg_ranking_row['stoic']

    print("\n!!!!!!!!!!Final Selection!!!!!!!!!!\n")
    print(f"Stoichiometry with the highest Maximum ranking score: {stoic_max_ranking}")
    print(f"Stoichiometry with the highest Average ranking score: {stoic_avg_ranking}")


if __name__ == '__main__':
    """
    --input_fasta : *.fasta file location 
    --stoichiometries : Comma separated stoichiometries. Example --stoichiometries A2,A3,A4
    --num_models : Number of models to generate for each stoichiometry (in the multiple of 5). Example --num_model 25  
    --output_path : Directory where output is to be saved. Example --output_path --/home/user/examples/
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta',type=is_file,required=True)
    parser.add_argument('--stoichiometries', required=True)
    parser.add_argument('--num_models', required=True)
    parser.add_argument('--output_path', type=is_dir, required=True)

    args = parser.parse_args()
    ## Get Alphafold3 configurations
    if not os.path.exists("config.json"):
        print("Please run configure.py first to create a config.json file containing af3_program_path, af3_params_path, af3_db_path, hhsearch_binary, hhmake_binary, and hhdb_prefix")
        exit()
    with open("config.json", "r") as f:
        config = json.load(f)
    af3_program_path = config["af3_program_path"]
    af3_params_path = config["af3_params_path"]
    af3_db_path = config["af3_db_path"]

    # Receive and validate input information
    target_name, sequences = read_fasta(args.input_fasta)
    stoichiometries = args.stoichiometries.split(",")
    stoichiometries = [i.lower() for i in stoichiometries]
    unique_chains = process_stoichiometry(stoichiometries)
    default_stoichiometry = create_default_stoichiometry(unique_chains)
    if default_stoichiometry not in stoichiometries: stoichiometries.append(default_stoichiometry)
    stoichiometries = sorted(stoichiometries)
    if not validate_sequence_and_stoichiometry(sequences,unique_chains):
        print("Error in the provided information")
        exit()
    num_seeds = int(int(args.num_models)/5)# each seed contributes to five models so

    #Setup output paths
    output_path = args.output_path
    makedir_if_not_exists(output_path)
    target_output_path = os.path.join(output_path,target_name)

    # Generate the required models for different stoichiometries using AlphaFold3 program.
    generate_structures(af3_program_path,af3_params_path,af3_db_path,target_name,target_output_path,sequences,stoichiometries,default_stoichiometry,num_seeds)
    
    # Print the results for the models generated for different stoichiometries above.
    print_result(target_name,target_output_path,stoichiometries)

    

    
