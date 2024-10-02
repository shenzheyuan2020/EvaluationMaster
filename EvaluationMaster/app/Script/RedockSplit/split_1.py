# import os
# import sys
# import numpy as np
# import subprocess
# import csv
# from Bio.PDB import PDBParser, PDBIO, Select
# from rdkit import Chem
# from rdkit.Chem import AllChem
# # import pdbfixer
# # from openmm.app import PDBFile
# import logging

# # 定义氨基酸残基名称
# standard_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
#                      'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
#                      'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
# water_resnames = ['HOH', 'WAT', 'H2O']
# metal_resnames = ['NA', 'K', 'MG', 'CA', 'MN', 'FE', 'ZN', 'CU', 'CO', 'NI', 'CD', 'HG', 'SO4']
# modified_resnames = ['PTR', 'SEP', 'TPO']

# def is_het(residue):
#     """判断是否为异质原子（非标准残基，如水和配体）"""
#     res = residue.id[0]
#     return res != " " and res != "W"

# def is_removable_residue(residue):
#     """判断是否为需要删除的残基，包括水分子、金属离子、酸根离子和修饰后的氨基酸"""
#     resname = residue.get_resname().strip()
#     return (resname in water_resnames or
#             resname in metal_resnames or
#             resname in modified_resnames)

# class ResidueSelect(Select):
#     def __init__(self, chain, residue):
#         self.chain = chain
#         self.residue = residue

#     def accept_chain(self, chain):
#         return chain.id == self.chain.id

#     def accept_residue(self, residue):
#         return residue == self.residue and is_het(residue)

# class AtomOnlySelect(Select):
#     """只选择蛋白质的ATOM行，不包含HETATM行（如配体和其他异质原子）"""
#     def accept_residue(self, residue):
#         return residue.get_resname() in standard_resnames

# def preprocess_structure(structure):
#     """删除结构中的水分子、金属离子、酸根离子和修饰后的氨基酸"""
#     for model in structure:
#         for chain in model:
#             residues_to_remove = []
#             for residue in chain:
#                 if is_removable_residue(residue):
#                     residues_to_remove.append(residue.id)

#             for res_id in residues_to_remove:
#                 chain.detach_child(res_id)

# def extract_unique_ligand(structure):
#     """提取唯一的配体及其对应的蛋白链"""
#     unique_ligand = None
#     unique_chain = None

#     ligands = []
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if is_het(residue) and residue.get_resname() not in water_resnames:
#                     ligands.append((chain, residue))

#     if ligands:
#         unique_chain, unique_ligand = ligands[0]

#     return unique_chain, unique_ligand

# def generate_smiles_from_residue(chain, residue):
#     """生成残基的SMILES字符串"""
#     io = PDBIO()
#     ligand_pdb = f"temp_lig_{chain.id}.pdb"
#     io.set_structure(chain)
#     io.save(ligand_pdb, ResidueSelect(chain, residue))
#     mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
#     os.remove(ligand_pdb)
#     if mol:
#         mol = Chem.AddHs(mol)
#         return Chem.MolToSmiles(mol, canonical=True)
#     return None

# def save_origin_pdb(structure, output_file):
#     """保存处理后的结构为{pdbname}_origin.pdb文件"""
#     io = PDBIO()
#     io.set_structure(structure)
#     io.save(output_file)

# def extract_protein_chain_and_ligand(path, structure, pdbname, unique_chain, unique_ligand):
#     """提取蛋白链和配体，生成单独的文件"""
#     io = PDBIO()
    
#     # 保存蛋白链文件
#     protein_chain_file = os.path.join(path, f"{pdbname}_chain.pdb")
#     io.set_structure(unique_chain)
#     io.save(protein_chain_file, AtomOnlySelect())

#     # 保存配体文件
#     ligand_pdb = os.path.join(path, f"lig_{unique_chain.id}.pdb")
#     io.save(ligand_pdb, ResidueSelect(unique_chain, unique_ligand))

#     # 生成MOL2文件
#     reference_mol2 = os.path.join(path, "reference.mol2")
#     smiles_csv = os.path.join(path, "ligands_smiles.csv")

#     # 确保CSV文件存在并写入表头
#     if not os.path.exists(smiles_csv):
#         with open(smiles_csv, 'w', newline='') as f:
#             writer = csv.writer(f)
#             writer.writerow(["smiles", "Name"])  # 写入列名

#     smiles_file = convert_and_optimize_ligand(ligand_pdb, reference_mol2)

#     # 将SMILES写入CSV文件
#     if smiles_file:
#         save_ligand_smiles_csv(smiles_csv, smiles_file)

#     # 删除配体的PDB文件
#     if os.path.exists(ligand_pdb):
#         os.remove(ligand_pdb)

#     return protein_chain_file, smiles_csv  # 返回CSV路径

# def calculate_ligand_center(residue):
#     """计算配体的中心坐标"""
#     coords = np.mean([atom.get_coord() for atom in residue], axis=0)
#     return coords

# def save_pocket_center_csv(output_dir, protein_file, ligand_coords):
#     """保存蛋白口袋的中心坐标CSV文件"""
#     protein_file_base = os.path.splitext(os.path.basename(protein_file))[0]  # 去掉.pdb后缀
#     center_csv = os.path.join(output_dir, "pocket_center.csv")
#     with open(center_csv, 'w', newline='') as f:
#         writer = csv.writer(f)
#         writer.writerow(["mol_name", "x", "y", "z"])
#         writer.writerow([protein_file_base, ligand_coords[0], ligand_coords[1], ligand_coords[2]])

# def save_ligand_smiles_csv(smiles_csv, smiles_file):
#     """保存配体的SMILES字符串和名称到CSV文件"""
#     with open(smiles_file, 'r') as smi_f:
#         smiles_lines = smi_f.readlines()
    
#     with open(smiles_csv, 'a', newline='') as f:
#         writer = csv.writer(f)
#         for line in smiles_lines:
#             smiles, name = line.strip().split()  # 假设.smi文件中每行是"SMILES NAME"
#             writer.writerow([smiles, name])  # 写入SMILES和名称

# def convert_and_optimize_ligand(ligand_pdb, output_mol2):
#     mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=True)
#     if mol:
#         if AllChem.MMFFHasAllMoleculeParams(mol):
#             AllChem.MMFFOptimizeMolecule(mol)

#         # 优化后的文件转为MOL2格式
#         command = f"obabel {ligand_pdb} -O {output_mol2} --gen3d -r --ff MMFF94"
#         subprocess.run(command, shell=True)

#         # 生成标准SMILES并写入.smi文件
#         smiles_file = output_mol2.replace('.mol2', '.smi')
#         command = f"obabel {output_mol2} -O {smiles_file} --canonical"
#         subprocess.run(command, shell=True)

#         return smiles_file  # 返回.smi文件路径

#     return None

# def prepare_ligands_from_smiles(smiles_csv, base_output_dir, mgltool_path):
#     """
#     调用生成3D结构的脚本，将SMILES文件转换为SDF、MOL2和PDBQT格式。
    
#     :param smiles_csv: SMILES的CSV文件路径
#     :param base_output_dir: 输出文件夹路径
#     :param mgltool_path: MGLTools的安装路径
#     """
#     evaluation_master = os.getenv('EVALUATIONMASTER', '') #Get environment variable
#     ligprep_script = evaluation_master + "/EvaluationMaster/app/Script/LigScript/Gen3D.py"
#     lig_start_name = "ligand"  # 配体名称前缀
#     pdbqt_output_dir = os.path.join(base_output_dir, 'ligands', "PDBQT")  # 定义PDBQT输出目录
#     mol2_output_dir = os.path.join(base_output_dir, 'ligands', "MOL2")
#     sdf_output_dir = os.path.join(base_output_dir, 'ligands', "SDF")

#     # 创建输出目录
#     os.makedirs(pdbqt_output_dir, exist_ok=True)
#     os.makedirs(mol2_output_dir, exist_ok=True)
#     os.makedirs(sdf_output_dir, exist_ok=True)

#     try:
#         subprocess.run(['python', ligprep_script, lig_start_name, smiles_csv, base_output_dir, mgltool_path, '16'], check=True)
#         logging.info(f"Ligands prepared successfully from {smiles_csv} into {pdbqt_output_dir}, {mol2_output_dir} and {sdf_output_dir}")
#     except subprocess.CalledProcessError as e:
#         logging.error(f"Error during ligand preparation: {e}")
#         sys.exit(1)

# def process_protein(pdb_file, output_dir, mgltool_path):
#     pdbname = os.path.splitext(os.path.basename(pdb_file))[0]
#     split_data_dir = os.path.join(output_dir, f'split_data_{pdbname}')
#     os.makedirs(split_data_dir, exist_ok=True)

#     evaluation_master = os.getenv('EVALUATIONMASTER', '')
#     mgltool_paths = os.path.abspath(os.path.join(evaluation_master, mgltool_path))

#     print(f"Parsing PDB file: {pdb_file}")
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure(pdbname, pdb_file)

#     print("Preprocessing structure: removing water, ions, and modified residues...")
#     preprocess_structure(structure)

#     # 提取唯一的配体和蛋白链
#     print("Extracting unique ligand and corresponding chain...")
#     unique_chain, unique_ligand = extract_unique_ligand(structure)

#     # 保存origin.pdb文件
#     origin_pdb_file = os.path.join(split_data_dir, f"{pdbname}_origin.pdb")
#     save_origin_pdb(structure, origin_pdb_file)
#     print(f"Saved origin PDB file as: {origin_pdb_file}")

#     # 提取蛋白链和配体，生成对应文件
#     print("Extracting protein chain and ligand files...")
#     protein_chain_file, smiles_csv = extract_protein_chain_and_ligand(split_data_dir, structure, pdbname, unique_chain, unique_ligand)

#     # 计算配体中心并保存口袋中心坐标CSV
#     ligand_coords = calculate_ligand_center(unique_ligand)
#     save_pocket_center_csv(split_data_dir, protein_chain_file, ligand_coords)
#     print(f"Saved pocket center CSV file.")

#     # 准备配体
#     prepare_ligands_from_smiles(smiles_csv, split_data_dir, mgltool_paths)

# def main():
#     if len(sys.argv) != 4:
#         print("Usage: python script.py protein_dir output_dir mgltool_path")
#         sys.exit(1)

#     protein_dir = sys.argv[1]
#     output_dir = sys.argv[2]
#     mgltool_path = sys.argv[3]

#     for pdb_file in os.listdir(protein_dir):
#         if pdb_file.endswith(".pdb"):
#             pdb_file_path = os.path.join(protein_dir, pdb_file)
#             process_protein(pdb_file_path, output_dir, mgltool_path)

# if __name__ == "__main__":
#     main()

import os
import sys
import numpy as np
import subprocess
import csv
from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

# 定义氨基酸残基名称
standard_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
                     'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                     'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
water_resnames = ['HOH', 'WAT', 'H2O']
metal_resnames = ['NA', 'K', 'MG', 'CA', 'MN', 'FE', 'ZN', 'CU', 'CO', 'NI', 'CD', 'HG', 'SO4']
modified_resnames = ['PTR', 'SEP', 'TPO']

def is_het(residue):
    """判断是否为异质原子（非标准残基，如水和配体）"""
    res = residue.id[0]
    return res != " " and res != "W"

def is_removable_residue(residue):
    """判断是否为需要删除的残基，包括水分子、金属离子、酸根离子和修饰后的氨基酸"""
    resname = residue.get_resname().strip()
    return (resname in water_resnames or
            resname in metal_resnames or
            resname in modified_resnames)

class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue

    def accept_chain(self, chain):
        return chain.id == self.chain.id

    def accept_residue(self, residue):
        return residue == self.residue and is_het(residue)

class AtomOnlySelect(Select):
    """只选择蛋白质的ATOM行，不包含HETATM行（如配体和其他异质原子）"""
    def accept_residue(self, residue):
        return residue.get_resname() in standard_resnames

def preprocess_structure(structure):
    """删除结构中的水分子、金属离子、酸根离子和修饰后的氨基酸"""
    for model in structure:
        for chain in model:
            residues_to_remove = []
            for residue in chain:
                if is_removable_residue(residue):
                    residues_to_remove.append(residue.id)

            for res_id in residues_to_remove:
                chain.detach_child(res_id)

def extract_unique_ligand(structure):
    """提取唯一的配体及其对应的蛋白链"""
    unique_ligand = None
    unique_chain = None

    ligands = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_het(residue) and residue.get_resname() not in water_resnames:
                    ligands.append((chain, residue))

    if ligands:
        unique_chain, unique_ligand = ligands[0]

    return unique_chain, unique_ligand

def save_origin_pdb(structure, output_file):
    """保存处理后的结构为{pdbname}_origin.pdb文件"""
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

def extract_protein_chain_and_ligand(path, structure, pdbname, unique_chain, unique_ligand):
    """提取蛋白链和配体，生成单独的文件"""
    io = PDBIO()
    
    # 保存蛋白链文件
    protein_chain_file = os.path.join(path, f"{pdbname}_chain.pdb")
    io.set_structure(unique_chain)
    io.save(protein_chain_file, AtomOnlySelect())

    # 保存配体文件
    ligand_pdb = os.path.join(path, f"{pdbname}_lig.pdb")
    io.save(ligand_pdb, ResidueSelect(unique_chain, unique_ligand))
    return ligand_pdb, protein_chain_file  # 返回配体和蛋白链的路径

def calculate_ligand_center(residue):
    """计算配体的中心坐标"""
    coords = np.mean([atom.get_coord() for atom in residue], axis=0)
    return coords

def save_pocket_center_csv(output_dir, protein_file, ligand_coords):
    """保存蛋白口袋的中心坐标CSV文件"""
    protein_file_base = os.path.splitext(os.path.basename(protein_file))[0]  # 去掉.pdb后缀
    center_csv = os.path.join(output_dir, "pocket_center.csv")
    with open(center_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["mol_name", "x", "y", "z"])
        writer.writerow([protein_file_base, ligand_coords[0], ligand_coords[1], ligand_coords[2]])

def convert_ligand_formats(ligand_pdb, output_dir, pdbname, mgltool_path):
    """使用RDKit将配体文件转换为SDF和SMILES格式，使用Open Babel转换为MOL2格式，再用MGLTools转换为PDBQT格式"""
    mol = Chem.MolFromPDBFile(ligand_pdb)
    if not mol:
        print(f"Error: Unable to read PDB file {ligand_pdb} using RDKit.")
        return

    # 生成SDF文件
    sdf_file = os.path.join(output_dir, f"{pdbname}_lig.sdf")
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()
    print(f"Saved SDF file: {sdf_file}")

    # 生成MOL2文件使用Open Babel
    mol2_file = os.path.join(output_dir, f"{pdbname}_lig.mol2")
    subprocess.run(['obabel', '-i', 'pdb', ligand_pdb, '-o', 'mol2', '-O', mol2_file])
    print(f"Saved MOL2 file using Open Babel: {mol2_file}")

    # 检查MOL2文件是否存在
    if not os.path.exists(mol2_file):
        print(f"Error: MOL2 file {mol2_file} does not exist.")
        return

    # 生成SMILES字符串
    smiles_file = os.path.join(output_dir, f"{pdbname}_lig.smiles")
    smiles = Chem.MolToSmiles(mol)
    with open(smiles_file, 'w') as f:
        f.write(smiles)
    print(f"Saved SMILES file: {smiles_file}")

    # 生成PDBQT文件使用MGLTools
    pdbqt_file = os.path.join(output_dir, f"{pdbname}_lig.pdbqt")
    convert_mol2_to_pdbqt(mol2_file, pdbqt_file, mgltool_path)
    print(f"Saved PDBQT file using MGLTools: {pdbqt_file}")

def convert_mol2_to_pdbqt(mol2_path, pdbqt_path, mgltool_path):
    """
    使用MGLTOOLS将mol2文件转换为PDBQT文件。

    :param mol2_path: MOL2文件的路径。
    :param pdbqt_path: 输出的PDBQT文件的路径。
    :param mgltool_path: MGLTOOLS的路径。
    """
    try:
        # 保存当前目录
        current_dir = os.getcwd()
        
        # 切换到MOL2文件所在目录
        mol2_dir = os.path.dirname(mol2_path)
        os.chdir(mol2_dir)
        
        # 执行转换
        subprocess.run([
            os.path.join(mgltool_path, 'bin', 'pythonsh'),
            os.path.join(mgltool_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_ligand4.py'),
            '-l', os.path.basename(mol2_path),
            '-o', os.path.basename(pdbqt_path),
            '-A', 'hydrogens'
        ], check=True)
        
        # 切换回原目录
        os.chdir(current_dir)
        
        logging.info(f"Successfully converted {mol2_path} to {pdbqt_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in converting {mol2_path} to PDBQT: {e}")
        # 切换回原目录，即使发生错误也要切换回来
        os.chdir(current_dir)

def process_protein(pdb_file, redock_dir, mgltool_path):
    pdbname = os.path.splitext(os.path.basename(pdb_file))[0]
    split_data_dir = os.path.join(redock_dir, f'split_data_{pdbname}')
    os.makedirs(split_data_dir, exist_ok=True)

    print(f"Parsing PDB file: {pdb_file}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdbname, pdb_file)

    print("Preprocessing structure: removing water, ions, and modified residues...")
    preprocess_structure(structure)

    # 提取唯一的配体和蛋白链
    print("Extracting unique ligand and corresponding chain...")
    unique_chain, unique_ligand = extract_unique_ligand(structure)

    # 保存origin.pdb文件
    origin_pdb_file = os.path.join(split_data_dir, f"{pdbname}_origin.pdb")
    save_origin_pdb(structure, origin_pdb_file)
    print(f"Saved origin PDB file as: {origin_pdb_file}")

    # 提取蛋白链和配体，生成对应文件
    print("Extracting protein chain and ligand files...")
    ligand_pdb, protein_chain_file = extract_protein_chain_and_ligand(split_data_dir, structure, pdbname, unique_chain, unique_ligand)

    # 计算配体中心并保存口袋中心坐标CSV
    ligand_coords = calculate_ligand_center(unique_ligand)
    save_pocket_center_csv(split_data_dir, protein_chain_file, ligand_coords)
    print(f"Saved pocket center CSV file.")

    # 转换配体格式
    print("Converting ligand to SDF, MOL2, SMILES, and PDBQT formats...")
    convert_ligand_formats(ligand_pdb, split_data_dir, pdbname, mgltool_path)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py protein_dir mgltool_path")
        sys.exit(1)

    protein_dir = sys.argv[1]
    mgltool_path = sys.argv[2]

    for pdb_file in os.listdir(protein_dir):
        if pdb_file.endswith(".pdb"):
            pdb_file_path = os.path.join(protein_dir, pdb_file)
            redock_dir = os.path.join(protein_dir, 'Redock')
            os.makedirs(redock_dir, exist_ok=True)
            process_protein(pdb_file_path, redock_dir, mgltool_path)

if __name__ == "__main__":
    main()
