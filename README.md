# prestoi
Predicting stoichiometry of protein complexes using AlphaFold3 and structural templates

# Begin with the installation of AlphaFold3 program using 
https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md

Test whether AlphaFold3 program is working properly
Once you have installed AlphaFold 3, you can test your setup using e.g. the
following input JSON file named `fold_input.json`:

```json
{
  "name": "2PV7",
  "sequences": [
    {
      "protein": {
        "id": ["A", "B"],
        "sequence": "GMRESYANENQFGFKTINSDIHKIVIVGGYGKLGGLFARYLRASGYPISILDREDWAVAESILANADVVIVSVPINLTLETIERLKPYLTENMLLADLTSVKREPLAKMLEVHTGAVLGLHPMFGADIASMAKQVVVRCDGRFPERYEWLLEQIQIWGAKIYQTNATEHDHNMTYIQALRHFSTFANGLHLSKQPINLANLLALSSPIYRLELAMIGRLFAQDAELYADIIMDKSENLAVIETLKQTYDEALTFFENNDRQGFIDAFHKVRDWFGDYSEQFLKESRQLLQQANDLKQG"
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
```

You can then run AlphaFold 3 using the following command:

```
docker run -it \
    --volume $HOME/af_input:/root/af_input \
    --volume $HOME/af_output:/root/af_output \
    --volume <MODEL_PARAMETERS_DIR>:/root/models \
    --volume <DATABASES_DIR>:/root/public_databases \
    --gpus all \
    alphafold3 \
    python run_alphafold.py \
    --json_path=/root/af_input/fold_input.json \
    --model_dir=/root/models \
    --output_dir=/root/af_output
```

# Steps to run the stoichiometry prediction:
- Copy the codes to alphafold3/ directory
  ```
  cp stoichiometry_prediction.py protein_utils.py utils.py /path/to/alphafold3
  cd /path/to/alphafold3
  ```
- Run the stoichiometry_prediction.py
## Homomultimer Example
```
python stoichiometry_prediction.py --input_fasta /path/to/input_fasta --stoichiometries A2,A3,A4 --output_path /path/to/output_dir --db_path /bmlfast/databases/ --params_path /path/to/alphafold3_databases --num_models 25
```
## Heteromultimer Example
```
python stoichiometry_prediction.py --input_fasta /path/to/input_fasta --stoichiometries A1B1,A2B2,A9B18 --output_path /path/to/output_dir --db_path /bmlfast/databases/ --params_path /path/to/alphafold3_databases --num_models 25
```



