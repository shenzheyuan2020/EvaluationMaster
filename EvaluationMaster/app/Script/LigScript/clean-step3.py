# import sys

# import pandas as pd

# from rdkit import Chem

# from rdkit.Chem import SaltRemover

# import numpy as np

# import matplotlib.pyplot as plt

# import os

# from tqdm import tqdm



# # def clean_smiles(smiles_list, max_atoms=80):

# #     clean_smiles = []

# #     remover = SaltRemover.SaltRemover()

# #     for smiles in tqdm(smiles_list):

# #         try:

# #             mol = Chem.MolFromSmiles(smiles)

# #             if mol is not None and mol.GetNumAtoms() <= max_atoms:

# #                 mol = remover.StripMol(mol)

# #                 clean_smiles.append(Chem.MolToSmiles(mol))

# #         except:

# #             continue

# #     return clean_smiles

# def clean_smiles(smiles_list, max_atoms=80):
#     clean_smiles_pairs = []  # 存储原始smiles和清洗后的smiles对
#     remover = SaltRemover.SaltRemover()
#     for smiles in tqdm(smiles_list):
#         try:
#             mol = Chem.MolFromSmiles(smiles)
#             if mol is not None and mol.GetNumAtoms() <= max_atoms:
#                 mol = remover.StripMol(mol)
#                 clean_smiles_pairs.append((smiles, Chem.MolToSmiles(mol)))  # 保存原始和清洗后的smiles
#         except:
#             continue
#     return clean_smiles_pairs


# def plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix=''):

#     plt.figure(figsize=(10, 8))

#     plt.hist(compounds_df[column].dropna(), bins=50, color='skyblue', edgecolor='slategray', alpha=0.7)

#     plt.grid(axis='y', alpha=0.75)

#     plt.xlabel(column)

#     plt.ylabel('Frequency')

#     title = f'{filter_name} Enhanced Histogram of {column} {title_suffix}'.strip()

#     plt.title(title)



#     pic50_to_ic50 = [

#         (5, 10000),

#         (6, 1000),

#         (7, 100),

#         (8, 10),

#         (9, 1) 

#     ]



#     textstr = '\n'.join([f'pIC50 = {pic50}: IC50 = {ic50} nM' for pic50, ic50 in pic50_to_ic50])

#     props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

#     plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,

#              verticalalignment='top', horizontalalignment='right', bbox=props)

    

#     image_file_path = os.path.join(save_dir, f"{filter_name}_clean.png")

#     plt.savefig(image_file_path, bbox_inches='tight')

#     plt.close()

#     print(f"The enhanced histogram and comparison table have been saved to: {image_file_path}")



# def main(csv_file_path, filter_name, save_dir):

#     compounds_df = pd.read_csv(csv_file_path)

#     cleaned_smiles = clean_smiles(compounds_df['smiles'])

#     cleaned_df = pd.DataFrame(cleaned_smiles, columns=['smiles'])

#     cleaned_df = cleaned_df.drop_duplicates()



#     compounds_df = compounds_df.iloc[:len(cleaned_df)].copy()

#     compounds_df['smiles'] = cleaned_df['smiles'].values

    

#     cleaned_csv_file_path = os.path.join(save_dir, f'{filter_name}_clean.csv')

#     compounds_df.to_csv(cleaned_csv_file_path, index=False)

#     print(f"Cleaned data saved to: {cleaned_csv_file_path}")



#     plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix='_cleaned')

# def main(csv_file_path, filter_name, save_dir):
#     compounds_df = pd.read_csv(csv_file_path)
#     clean_smiles_pairs = clean_smiles(compounds_df['smiles'])  # 调整后的函数返回值
    
#     cleaned_df = pd.DataFrame(clean_smiles_pairs, columns=['original_smiles', 'cleaned_smiles'])
#     cleaned_df = cleaned_df.drop_duplicates(subset='cleaned_smiles')
    
#     # 更新compounds_df以包含清洗后的smiles
#     compounds_df['smiles'] = [pair[1] for pair in clean_smiles_pairs][:len(compounds_df)]
    
#     # 保存清洗后的数据为CSV
#     cleaned_csv_file_path = os.path.join(save_dir, f'{filter_name}_clean.csv')
#     compounds_df.to_csv(cleaned_csv_file_path, index=False)
#     print(f"Cleaned data saved to: {cleaned_csv_file_path}")
    
#     # 保存为.smi文件
#     smi_file_path = os.path.join(save_dir, f'{filter_name}_clean.smi')
#     with open(smi_file_path, 'w') as smi_file:
#         for _, row in compounds_df.iterrows():
#             smi_file.write(f"{row['cleaned_smiles']} {row['chembl_id']}\n")
#     print(f"SMI file has been saved to: {smi_file_path}")
    
#     plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix='_cleaned')




# if __name__ == "__main__":

#     if len(sys.argv) != 4:

#         print("Usage: python script.py <csv_file_path> <filter_name> <save_dir>")

#         sys.exit(1)



#     csv_file_path = sys.argv[1]

#     filter_name = sys.argv[2]

#     save_dir = sys.argv[3]

#     main(csv_file_path, filter_name, save_dir)

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SaltRemover
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

def clean_smiles(smiles_list, max_atoms=80):
    clean_smiles_pairs = []  # 存储原始smiles和清洗后的smiles对
    remover = SaltRemover.SaltRemover()
    for smiles in tqdm(smiles_list):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None and mol.GetNumAtoms() <= max_atoms:
                mol = remover.StripMol(mol)
                clean_smiles_pairs.append(Chem.MolToSmiles(mol))  # 只保存清洗后的smiles
            else:
                clean_smiles_pairs.append(None)  # 如果不符合条件，添加 None
        except:
            clean_smiles_pairs.append(None)  # 处理异常时也添加 None
    return clean_smiles_pairs

def plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix=''):
    plt.figure(figsize=(10, 8))
    plt.hist(compounds_df[column].dropna(), bins=50, color='skyblue', edgecolor='slategray', alpha=0.7)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(column)
    plt.ylabel('Frequency')
    title = f'{filter_name} Enhanced Histogram of {column} {title_suffix}'.strip()
    plt.title(title)
    
    # 根据不同的列动态生成对照表
    if column == 'pIC50':
        pic50_to_ic50 = [
            (5, 10000),
            (6, 1000),
            (7, 100),
            (8, 10),
            (9, 1)  
        ]
        textstr = '\n'.join([f'pIC50 = {pic50}: IC50 = {ic50} nM' for pic50, ic50 in pic50_to_ic50])
    elif column == 'pKd':
        kd_to_pkd = [
            (5, 10000),
            (6, 1000),
            (7, 100),
            (8, 10),
            (9, 1)
        ]
        textstr = '\n'.join([f'pKd = {pkd}: Kd = {kd} nM' for pkd, kd in kd_to_pkd])
    elif column == 'pKi':
        ki_to_pki = [
            (5, 10000),
            (6, 1000),
            (7, 100),
            (8, 10),
            (9, 1)
        ]
        textstr = '\n'.join([f'pKi = {pki}: Ki = {ki} nM' for pki, ki in ki_to_pki])

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right', bbox=props)
    
    image_file_path = os.path.join(save_dir, f"{filter_name}_clean.png")
    plt.savefig(image_file_path, bbox_inches='tight')
    plt.close()
    print(f"The {column} histogram has been saved to: {image_file_path}")

def main(csv_file_path, filter_name, save_dir):
    compounds_df = pd.read_csv(csv_file_path)

    # 清洗 SMILES 数据并返回清洗后的结果
    clean_smiles_list = clean_smiles(compounds_df['smiles'])
    compounds_df['smiles'] = clean_smiles_list

    # 移除清洗失败的行
    compounds_df = compounds_df[compounds_df['smiles'].notna()]
    
    # 保存清洗后的数据为 CSV
    cleaned_csv_file_path = os.path.join(save_dir, f'{filter_name}_clean.csv')
    compounds_df.to_csv(cleaned_csv_file_path, index=False)
    print(f"Cleaned data saved to: {cleaned_csv_file_path}")
    
    # 保存为 SMI 文件
    smi_file_path = os.path.join(save_dir, f'{filter_name}_clean.smi')
    with open(smi_file_path, 'w') as smi_file:
        for _, row in compounds_df.iterrows():
            smi_file.write(f"{row['cleaned_smiles']} {row.get('chembl_id', '')}\n")
    print(f"SMI file has been saved to: {smi_file_path}")

    # 计算 pIC50、pKd 和 pKi
    if 'ic50' in compounds_df.columns:
        compounds_df['pIC50'] = -np.log10(compounds_df['ic50'] * 1e-9)
    if 'kd' in compounds_df.columns:
        compounds_df['pKd'] = -np.log10(compounds_df['kd'] * 1e-9)
    if 'ki' in compounds_df.columns:
        compounds_df['pKi'] = -np.log10(compounds_df['ki'] * 1e-9)

    # 画图
    if 'pIC50' in compounds_df.columns:
        plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix='_cleaned')
    if 'pKd' in compounds_df.columns:
        plot_distribution(compounds_df, filter_name, save_dir, column='pKd', title_suffix='_cleaned')
    if 'pKi' in compounds_df.columns:
        plot_distribution(compounds_df, filter_name, save_dir, column='pKi', title_suffix='_cleaned')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <csv_file_path> <filter_name> <save_dir>")
        sys.exit(1)

    csv_file_path = sys.argv[1]
    filter_name = sys.argv[2]
    save_dir = sys.argv[3]
    main(csv_file_path, filter_name, save_dir)
