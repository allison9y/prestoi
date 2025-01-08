import os
import argparse
import json

def is_dir(path):
    if os.path.isdir(path):
        return os.path.abspath(path)
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid directory.")
def main():
    parser = argparse.ArgumentParser(description="Generate a JSON configuration file for Alphafold3 program.")
    parser.add_argument('--af3_program_path', type=is_dir, required=True, help="Path to the AlphaFold3 program directory.")
    parser.add_argument('--af3_params_path', type=is_dir, required=True, help="Path to the AlphaFold3 model parameters directory.")
    parser.add_argument('--af3_db_path', type=is_dir, required=True, help="Path to the AlphaFold3 database directory.")
    args = parser.parse_args()

    config = {
        "af3_program_path": args.af3_program_path,
        "af3_params_path": args.af3_params_path,
        "af3_db_path": args.af3_db_path
    }

    config_file = "config.json"
    with open(config_file, "w") as f:
        json.dump(config, f, indent=4)

    print(f"Configuration file '{config_file}' created successfully!")

if __name__ == "__main__":
    main()
