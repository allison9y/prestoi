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


def run_docker( op_path, params_path, db_path, target_name):
    docker_command = [
        "docker", "run", "--rm",  # Added --rm to automatically clean up container after execution
        "--volume", f"{op_path}:/root/af_output",
        "--volume", f"{params_path}:/root/models",
        "--volume", f"{db_path}:/root/public_databases",
        "--gpus", "all",
        "alphafold3",
        "python", "run_alphafold.py",
        f"--json_path=/root/af_output/input_jsons/{target_name}.json",
        "--model_dir=/root/models",
        "--output_dir=/root/af_output"
    ]
    print("Running Docker Command : ", " ".join(docker_command))

    try:
        process = subprocess.Popen(
            docker_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )

        for stdout_line in iter(process.stdout.readline, ""):
            print(stdout_line, end="")
        for stderr_line in iter(process.stderr.readline, ""):
            print(stderr_line, end="")

        process.stdout.close()
        process.stderr.close()
        return_code = process.wait()

        if return_code != 0:
            print(f"Docker command failed with return code {return_code}")
        else:
            print("Docker command completed successfully.")

    except Exception as e:
        print(f"Error occurred: {str(e)}")

    
def generate_structures(op_path,target_name,sequences,stoichiometries,default_stoichiometry,num_seeds):
    makedir_if_not_exists(op_path)
    ip_json_path = os.path.join(op_path,"input_jsons")
    makedir_if_not_exists(ip_json_path)

    generate_json(target_name,default_stoichiometry,num_seeds,ip_json_path,sequences,common_data={},default_stoichiometry=True)

    print(f"####################### Running for default Stoichiometry: {default_stoichiometry} #######################\n")
    
    if not os.path.exists(f"{op_path}/{target_name}_{default_stoichiometry}/ranking_scores.csv"):
        run_docker(op_path,params_path,db_path,f"{target_name}_{default_stoichiometry}")
    
    print(f"####################### Completed for Stoichiometry: {default_stoichiometry} #######################\n")

    common_data_path = f"{op_path}/{target_name}_{default_stoichiometry}/{target_name}_{default_stoichiometry}_data.json"
    if not os.path.exists(common_data_path):
        print("program Completed without generating results")
        print(common_data_path)
        exit()

    # Processing for other stoichiometries A3,A4......
    print(f"####################### Creating input.json for other stoichiometries #######################\n")

    common_data_dict = read_common_data(common_data_path)
    
    for stoic in stoichiometries:
        if stoic != default_stoichiometry:
            generate_json(target_name,stoic,num_seeds,ip_json_path,sequences,common_data=common_data_dict,default_stoichiometry=False)
        
    for stoic in stoichiometries[1:]:
        if stoic != default_stoichiometry:
            print(f"####################### Running for Stoichiometry: {stoic} #######################\n")
            ranking_csv_path = f"{op_path}/{target_name}_{str(stoic)}/ranking_scores.csv"
            if not os.path.exists(ranking_csv_path):
                run_docker(op_path,params_path,db_path,f"{target_name}_{stoic}")
            print(f"####################### Completed for Stoichiometry: {stoic} #######################\n")


def result_print(op_path,target_name,stoichiometries):
    print("Stoichiometry, Maximum ranking score, Average ranking score, Number of models") 
    data = [["stoic","max_ranking_score","avg_ranking_score","num_models"]]   
    for stoic in stoichiometries:
        ranking_csv_path = f"{op_path}/{target_name}_{str(stoic)}/ranking_scores.csv"
        if os.path.exists(ranking_csv_path):
            df = pd.read_csv(ranking_csv_path)
            ranking_score_list = df["ranking_score"].tolist()
            max_ranking_score = max(ranking_score_list)
            avg_ranking_score = sum(ranking_score_list)/len(ranking_score_list)
            model_count = len(ranking_score_list)
            print(f"{stoic.upper()},{max_ranking_score},{avg_ranking_score},{model_count}")
            data.append([stoic.upper(),max_ranking_score,avg_ranking_score,model_count])

    result_csv_path = os.path.join(op_path,'stoichiometry_results.csv')
    with open(result_csv_path,'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)

    df = pd.read_csv(result_csv_path)

    max_ranking_row = df.loc[df['max_ranking_score'].idxmax()]
    avg_ranking_row = df.loc[df['avg_ranking_score'].idxmax()]
    stoic_max_ranking = max_ranking_row['stoic']
    stoic_avg_ranking = avg_ranking_row['stoic']

    print("\n!!!!!!!!!!Final Results!!!!!!!!!!\n")
    print(f"Stoichiometry with highest Maximum ranking score: {stoic_max_ranking}")
    print(f"Stoichiometry with highest Average ranking score: {stoic_avg_ranking}")


if __name__ == '__main__':
    """
    --input_fasta : *.fasta file location 
    --stoichiometries : comma separated stoichiometries. Example --stoichiometries A2,A3,A4
    --num_models : number of models to generate for each stoichiometry (in the multiple of 5). Example --num_model 25  
    --output_path : directory where output is to be saved. Example --output_path --/bmlfast/pngkg/examples/
    --db_path : path to Alphafold3 databases. Example --db_path /bmlfast/databases
    --params_path : path to AlphaFold3 model parameters. Example --paramas_path /path/to/params
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta',type=is_file,required=True)
    parser.add_argument('--stoichiometries', required=True)
    parser.add_argument('--num_models', required=True)
    parser.add_argument('--output_path', type=is_dir, required=True)
    parser.add_argument('--db_path', type=is_dir, required=True)
    parser.add_argument('--params_path', type=is_dir, required=True)

    args = parser.parse_args()

    target_name, sequences = read_fasta(args.input_fasta)
    stoichiometries = args.stoichiometries.split(",")
    stoichiometries = [i.lower() for i in stoichiometries]
    unique_chains = process_stoichiometry(stoichiometries)
    default_stoichiometry = create_default_stoichiometry(unique_chains)
    if default_stoichiometry not in stoichiometries: stoichiometries.append(default_stoichiometry)
    stoichiometries = sorted(stoichiometries)
    if not validate_sequence_and_stoichiometry(sequences,unique_chains):
        exit()
    num_seeds = int(int(args.num_models)/5)# each seed contributes to five models so
    output_path = args.output_path
    makedir_if_not_exists(output_path)
    op_path = os.path.join(output_path,target_name)
    db_path = args.db_path
    params_path = args.params_path

    generate_structures(op_path,target_name,sequences,stoichiometries,default_stoichiometry,num_seeds)
    result_print(op_path,target_name,stoichiometries)

    

    