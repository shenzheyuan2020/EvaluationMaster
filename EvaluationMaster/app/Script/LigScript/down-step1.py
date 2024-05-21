# pip install pandas chembl_webresource_client matplotlib rdkit tqdm numpy
# python down.py uniport_ID name save_dir
import os
import sys
import pandas as pd
from tqdm.auto import tqdm
from chembl_webresource_client.new_client import new_client
import matplotlib.pyplot as plt
import numpy as np

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
    compounds_df = compounds_df.drop_duplicates(subset='smiles')
    return compounds_df

def convert_ic50_to_pic50(compounds_df):
    compounds_df['ic50'] = pd.to_numeric(compounds_df['ic50'], errors='coerce')
    compounds_df = compounds_df[compounds_df['ic50'] > 0]
    compounds_df.loc[:, 'pIC50'] = -np.log10(compounds_df['ic50'] * 1e-9)
    return compounds_df


def plot_distribution(compounds_df, name, save_dir, column='pIC50'):
    plt.figure(figsize=(10, 8))
    plt.hist(compounds_df[column].dropna(), bins=50, color='skyblue', edgecolor='slategray', alpha=0.7)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(column)
    plt.ylabel('Frequency')
    plt.title(f'{name} Enhanced Histogram of {column}')
    
    # 在直方图旁边添加对照表
    # 定义pIC50和IC50的对照列表
    pic50_to_ic50 = [
        (5, 10000),
        (6, 1000),
        (7, 100),
        (8, 10),
        (9, 1)  # pIC50值对应IC50值（单位：nM）
    ]
    
    # 添加文本框以显示对照表信息
    textstr = '\n'.join([f'pIC50 = {pic50}: IC50 = {ic50} nM' for pic50, ic50 in pic50_to_ic50])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right', bbox=props)
    
    image_file_path = os.path.join(save_dir, f"{name}_downlig.png")
    plt.savefig(image_file_path, bbox_inches='tight')
    plt.close()
    print(f"The enhanced histogram and comparison table have been saved to: {image_file_path}")



def main(uniprot_id, name, save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    compounds_df = get_active_compounds_from_chembl(uniprot_id)
    if compounds_df is not None:
        compounds_df = convert_ic50_to_pic50(compounds_df)
        csv_file_path = os.path.join(save_dir, f'{name}.csv')
        compounds_df.to_csv(csv_file_path, index=False)
        print(f"CSV file has been saved to: {csv_file_path}")
        plot_distribution(compounds_df, name, save_dir, column='pIC50')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <uniprot_id> <name> <save_dir>")
        sys.exit(1)

    uniprot_id = sys.argv[1]
    name = sys.argv[2]
    save_dir = sys.argv[3]
    main(uniprot_id, name, save_dir) 