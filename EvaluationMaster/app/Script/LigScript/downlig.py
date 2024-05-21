# pip install pandas chembl_webresource_client matplotlib rdkit tqdm numpy
import os
import sys
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import SaltRemover
import matplotlib.pyplot as plt


# 创建或确认 working_place/ 子目录的存在
cwd = os.getcwd()
# working_dir = os.path.join(cwd, 'working_place/lig/')
# if not os.path.exists(working_dir):
#     os.makedirs(working_dir)

def get_active_compounds_from_chembl(uniprot_id):
    if pd.isnull(uniprot_id):
        print(f"Uniprot ID is nan, skipping.")
        return None
    target = new_client.target
    target_query = target.search(uniprot_id)
    if target_query:
        target_id = target_query[0]['target_chembl_id']
    else:
        print(f"No target found for Uniprot ID: {uniprot_id}")
        return None
    activity = new_client.activity
    activity_query = activity.filter(target_chembl_id=target_id).filter(standard_type__in=["IC50"])

    compounds = []
    for act in tqdm(activity_query):
        compound = {
            'chembl_id': act['molecule_chembl_id'],
            'smiles': act['canonical_smiles'],
        }
        compound[act['standard_type'].lower()] = act['standard_value']
        compound[act['standard_type'].lower() + '_unit'] = act['standard_units'] 
        compounds.append(compound)
    compounds_df = pd.DataFrame(compounds)
    compounds_df = compounds_df.fillna(0)
    return compounds_df

def clean_smiles(smiles_list, max_atoms=80):
    clean_smiles = []
    remover = SaltRemover.SaltRemover()
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None and mol.GetNumAtoms() <= max_atoms:
            mol = remover.StripMol(mol)
            clean_smiles.append(Chem.MolToSmiles(mol))
    return clean_smiles

def process_uniprot_id(uniprot_id, name, save_dir):
    active_compounds = get_active_compounds_from_chembl(uniprot_id)
    if active_compounds is not None:
        smiles_list = active_compounds['smiles'].tolist()
        clean_smiles_list = clean_smiles(smiles_list)
        clean_compounds = active_compounds[active_compounds['smiles'].isin(clean_smiles_list)][['chembl_id', 'smiles', 'ic50', 'ic50_unit']]
        clean_compounds = clean_compounds.drop_duplicates(subset=['smiles'])
        clean_compounds['ic50'] = pd.to_numeric(clean_compounds['ic50'], errors='coerce')
        clean_compounds = clean_compounds.dropna(subset=['ic50'])
        clean_compounds = clean_compounds.sort_values(by='ic50')

        clean_compounds['ic50_M'] = clean_compounds['ic50'].astype(float) / 1_000_000_000
        clean_compounds = clean_compounds[clean_compounds['ic50_M'] != 0]
        clean_compounds['pIC50'] = -np.log10(clean_compounds['ic50_M'])
        clean_compounds = clean_compounds.drop(columns=['ic50_M'])

        csv_file_path = os.path.join(save_dir, f'{name}.csv')
        clean_compounds.to_csv(csv_file_path, index=False)
        print(f"CSV file has been saved to: {csv_file_path}")

        plt.figure(figsize=(10, 5))
        plt.figure(figsize=(10, 5))
        plt.hist(clean_compounds['pIC50'], bins=50, color='blue', alpha=0.7)
        plt.xlabel('pIC50')
        plt.ylabel('Number of molecules')
        plt.title(f'The distribution of pIC50 values for {name}')
        image_file_path = os.path.join(save_dir, f"{name}_pIC50_distribution.png")
        plt.savefig(image_file_path)
        print(f"The pIC50 distribution chart has been saved to: {image_file_path}")
        return csv_file_path, image_file_path
    return None, None

def main(uniprot_id, name, save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    process_uniprot_id(uniprot_id, name, save_dir)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <uniprot_id> <name> <save_dir>")
        sys.exit(1)

    uniprot_id = sys.argv[1]
    name = sys.argv[2]
    save_dir = sys.argv[3]
    main(uniprot_id, name, save_dir)
