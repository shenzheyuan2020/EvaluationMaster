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
    return center_csv  # 返回CSV文件的路径

def convert_ligand_formats(ligand_pdb, output_dir, pdbname, mgltool_path):
    """使用RDKit将配体文件转换为SDF和SMILES格式，使用Open Babel转换为MOL2格式，再用MGLTools转换为PDBQT格式"""
    mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
    if not mol:
        print(f"Error: Unable to read PDB file {ligand_pdb} using RDKit.")
        return None, None, None, None

    # 添加氢原子
    mol = Chem.AddHs(mol)

    # 生成SDF文件
    sdf_file = os.path.join(output_dir, f"{pdbname}_lig.sdf")
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()
    print(f"Saved SDF file: {sdf_file}")

    # 生成MOL2文件使用Open Babel
    mol2_file = os.path.join(output_dir, f"{pdbname}_lig.mol2")
    subprocess.run(['obabel', '-i', 'pdb', ligand_pdb, '-o', 'mol2', '-O', mol2_file, '-h'])
    print(f"Saved MOL2 file using Open Babel: {mol2_file}")

    # 检查MOL2文件是否存在
    if not os.path.exists(mol2_file):
        print(f"Error: MOL2 file {mol2_file} does not exist.")
        return None, None, None, None

    # 生成SMILES字符串并保存为CSV
    smiles_file = os.path.join(output_dir, f"{pdbname}_lig.smiles")
    smiles = Chem.MolToSmiles(mol)
    with open(smiles_file, 'w') as f:
        f.write(smiles)
    print(f"Saved SMILES file: {smiles_file}")

    smiles_csv = os.path.join(output_dir, f"{pdbname}_lig_smiles.csv")
    with open(smiles_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["smiles", "Name"])
        writer.writerow([smiles, f"{pdbname}_lig"])
    print(f"Saved SMILES CSV file: {smiles_csv}")

    # 生成PDBQT文件使用MGLTools
    pdbqt_file = os.path.join(output_dir, f"{pdbname}_lig.pdbqt")
    convert_mol2_to_pdbqt(mol2_file, pdbqt_file, mgltool_path)
    print(f"Saved PDBQT file using MGLTools: {pdbqt_file}")

    return sdf_file, mol2_file, smiles_csv, pdbqt_file  # 返回生成的文件路径

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

def calculate_rmsd_karmadock(karmadock_output_dir, original_ligand_file, pdbname):
    """
    计算 KarmaDock 输出的三个文件与原始配体之间的 RMSD，并将结果保存到 CSV 文件。

    :param karmadock_output_dir: KarmaDock 输出文件夹路径。
    :param original_ligand_file: 原始配体文件路径（PDB 格式）。
    :param pdbname: 蛋白名称，用于构建文件名。
    """
    # Paths to the KarmaDock output files
    files = [
        f"{pdbname}_lig_pred_align_corrected.sdf",
        f"{pdbname}_lig_pred_ff_corrected.sdf",
        f"{pdbname}_lig_pred_random_pose.sdf"
    ]
    methods = ["align_corrected", "ff_corrected", "random_pose"]
    rmsd_results = []

    for filename, method in zip(files, methods):
        sdf_path = os.path.join(karmadock_output_dir, filename)
        if os.path.exists(sdf_path):
            # Read the molecules
            mol_ref = Chem.MolFromPDBFile(original_ligand_file, removeHs=False)
            supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
            mol = supplier[0] if len(supplier) > 0 else None
            if mol_ref is None or mol is None:
                print(f"Error reading molecules for RMSD calculation: {sdf_path}")
                continue
            # Align the molecules and calculate RMSD
            rmsd = AllChem.GetBestRMS(mol_ref, mol)
            rmsd_results.append((method, rmsd))
        else:
            print(f"File not found: {sdf_path}")
    # Write the RMSD results to CSV
    rmsd_csv = os.path.join(karmadock_output_dir, f"{pdbname}_karmadock_rmsd.csv")
    with open(rmsd_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Method", "RMSD"])
        for method, rmsd in rmsd_results:
            writer.writerow([method, rmsd])
    print(f"Saved KarmaDock RMSD results to {rmsd_csv}")

def perform_docking(protein_pdb, center_csv, protein_pdb_dir, base_output_dir, pdbqt_output_dir, mol2_output_dir, smiles_csv, mgltool_path, autodock_gpu_path, karmadock_path, ledock_path, ligand_pdb, pdbname):
    """
    调用不同的对接工具（Vina, AutoDock-GPU, KarmaDock, LeDock）对蛋白和配体进行对接。

    :param protein_pdb: 处理后的蛋白质PDB文件路径
    :param center_csv: 配体中心坐标CSV文件路径
    :param protein_pdb_dir: 存放生成蛋白PDB的目录路径
    :param base_output_dir: 输出文件夹路径
    :param pdbqt_output_dir: PDBQT文件生成的目录路径
    :param mol2_output_dir: MOL2文件生成的目录路径
    :param smiles_csv: SMILES文件路径
    :param mgltool_path: MGLTools的路径
    :param autodock_gpu_path: AutoDock-GPU的路径
    :param karmadock_path: KarmaDock的路径
    :param ledock_path: LeDock的路径
    :param ligand_pdb: 原始配体PDB文件路径
    :param pdbname: 蛋白名称，用于构建文件名
    """
    logging.info(f"Starting docking process for PDB: {protein_pdb}")

    evaluation_master = os.getenv('EVALUATIONMASTER', '')

    # AutoDock Vina docking
    vina_gpu_script = f'{evaluation_master}/EvaluationMaster/app/Script/DockingScript/AutoDock_Vina.py'  # 替换为Vina对接脚本的路径
    try:
        subprocess.run([
            'python', vina_gpu_script, center_csv, pdbqt_output_dir, mgltool_path, base_output_dir
        ], check=True)
        logging.info(f"Vina-GPU docking completed for {protein_pdb}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in Vina-GPU docking: {e}")




    # Split Vina results
    dirnow = os.getcwd()
    pro_name = f'{os.path.splitext(os.path.basename(protein_pdb))[0]}_Vina'
    vina_output_dir = os.path.join(base_output_dir, pro_name)
    os.chdir(vina_output_dir)
    
    print(base_output_dir)
    vina_pdbqt = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig.pdbqt"
    print(f"Vina PDBQT file: {vina_pdbqt}")
    if os.path.exists(vina_pdbqt):
        try:
            subprocess.run(['vina_split', '--input', vina_pdbqt], check=True)
            logging.info(f"Vina results split for {vina_pdbqt}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in splitting Vina results: {e}")
    vina_pdbqtresult = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig_ligand_1.pdbqt"
    vina_convert = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig.sdf"
    subprocess.run(['obabel','-ipdbqt', vina_pdbqtresult, '-osdf', '-O', vina_convert])
    os.chdir(dirnow)
    # AutoDock-GPU docking
    autodock_gpu_script = f'{evaluation_master}/EvaluationMaster/app/Script/DockingScript/AutodockGPU.py'
    try:
        subprocess.run([
            'python', autodock_gpu_script, center_csv, pdbqt_output_dir, mgltool_path, autodock_gpu_path, base_output_dir, '2', '10'
        ], check=True)
        logging.info(f"AutoDock-GPU docking completed for {pdbqt_output_dir}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in AutoDock-GPU docking: {e}")


    dirnow = os.getcwd()
    pro_name = f'{os.path.splitext(os.path.basename(protein_pdb))[0]}_AutodockGPU'
    adgpu_output_dir = os.path.join(base_output_dir, pro_name)
    os.chdir(adgpu_output_dir)



    adgpu_pdbresult = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig-best.pdbqt"
    adgpu_convert = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig.sdf"
    subprocess.run(['obabel','-ipdbqt', adgpu_pdbresult, '-osdf', '-O', adgpu_convert])
    os.chdir(dirnow)

    # KarmaDock docking
    karmadock_script = f'{evaluation_master}/EvaluationMaster/app/Script/DockingScript/Karmadock.py'  # 替换为KarmaDock对接脚本路径
    try:
        subprocess.run([
            f'{karmadock_path}/Env/bin/python', karmadock_script, karmadock_path, center_csv, protein_pdb_dir, smiles_csv, base_output_dir, '0'
        ], check=True)
        logging.info(f"KarmaDock docking completed for {protein_pdb}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in KarmaDock docking: {e}")

    # Calculate RMSD for KarmaDock outputs
    pro_name = f'{os.path.splitext(os.path.basename(protein_pdb))[0]}_KarmaDock'
    karmadock_output_dir = os.path.join(base_output_dir, pro_name)
    calculate_rmsd_karmadock(karmadock_output_dir, ligand_pdb, pdbname)

    # LeDock docking
    ledock_script = f'{evaluation_master}/EvaluationMaster/app/Script/DockingScript/Ledock.py'  # 替换为LeDock对接脚本路径
    
    try:
        subprocess.run([
            'python', ledock_script, '4', center_csv, mol2_output_dir, base_output_dir, ledock_path
        ], check=True)
        logging.info(f"LeDock docking completed for {protein_pdb}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in LeDock docking: {e}")
    # Split LeDock results
    dirnow2 = os.getcwd()
    pro_name = f'{os.path.splitext(os.path.basename(protein_pdb))[0]}_ledock'
    ledock_output_dir = os.path.join(base_output_dir, pro_name)
    os.chdir(ledock_output_dir)


    ledock_dok = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig.dok"

    print(f"LeDock DOK file: {ledock_dok}")
    if os.path.exists(ledock_dok):
        try:
            subprocess.run([
                os.path.join(ledock_path, 'ledock'), '-spli', ledock_dok
            ], check=True)
            logging.info(f"LeDock results split for {ledock_dok}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in splitting LeDock results: {e}")

    ledock_pdbresult = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig_dock001.pdb"
    ledock_convert = f"{os.path.splitext(os.path.basename(protein_pdb))[0].rsplit('_', 1)[0]}_lig.sdf"
    subprocess.run(['obabel','-ipdb', ledock_pdbresult, '-osdf', '-O', ledock_convert])



    os.chdir(dirnow2)

def process_protein(pdb_file, redock_dir, mgltool_path, autodock_gpu_path, karmadock_path, ledock_path):
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

    if unique_chain is None or unique_ligand is None:
        print(f"No ligand found in {pdb_file}. Skipping.")
        return

    # 保存origin.pdb文件
    origin_pdb_file = os.path.join(split_data_dir, f"{pdbname}_origin.pdb")
    save_origin_pdb(structure, origin_pdb_file)
    print(f"Saved origin PDB file as: {origin_pdb_file}")

    # 提取蛋白链和配体，生成对应文件
    print("Extracting protein chain and ligand files...")
    ligand_pdb, protein_chain_file = extract_protein_chain_and_ligand(split_data_dir, structure, pdbname, unique_chain, unique_ligand)

    # 计算配体中心并保存口袋中心坐标CSV
    ligand_coords = calculate_ligand_center(unique_ligand)
    center_csv = save_pocket_center_csv(split_data_dir, protein_chain_file, ligand_coords)
    print(f"Saved pocket center CSV file: {center_csv}")

    # 转换配体格式
    print("Converting ligand to SDF, MOL2, SMILES, and PDBQT formats...")
    sdf_file, mol2_file, smiles_csv, pdbqt_file = convert_ligand_formats(ligand_pdb, split_data_dir, pdbname, mgltool_path)

    if None in (sdf_file, mol2_file, smiles_csv, pdbqt_file):
        print("Error in converting ligand formats. Skipping docking for this protein.")
        return

    # 收集需要的路径和参数
    protein_pdb = protein_chain_file
    protein_pdb_dir = os.path.dirname(protein_pdb)
    base_output_dir = split_data_dir
    pdbqt_output_dir = split_data_dir  # Assuming PDBQT files are in the same directory
    mol2_output_dir = split_data_dir   # Assuming MOL2 files are in the same directory

    # 调用对接函数
    perform_docking(
        protein_pdb=protein_pdb,
        center_csv=center_csv,
        protein_pdb_dir=protein_pdb_dir,
        base_output_dir=base_output_dir,
        pdbqt_output_dir=pdbqt_output_dir,
        mol2_output_dir=mol2_output_dir,
        smiles_csv=smiles_csv,
        mgltool_path=mgltool_path,
        autodock_gpu_path=autodock_gpu_path,
        karmadock_path=karmadock_path,
        ledock_path=ledock_path,
        ligand_pdb=ligand_pdb,
        pdbname=pdbname
    )


def main():
    if len(sys.argv) != 6:
        print("Usage: python script.py protein_dir mgltool_path autodock_gpu_path karmadock_path ledock_path")
        sys.exit(1)

    protein_dir = sys.argv[1]
    mgltool_path = sys.argv[2]
    autodock_gpu_path = sys.argv[3]
    karmadock_path = sys.argv[4]
    ledock_path = sys.argv[5]

    # Set environment variables if necessary
    evaluation_master = os.getenv('EVALUATIONMASTER', '')
    if not evaluation_master:
        print("Please set the EVALUATIONMASTER environment variable.")
        sys.exit(1)

    redock_dir = os.path.join(protein_dir, 'Redock')
    os.makedirs(redock_dir, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s')

    for pdb_file in os.listdir(protein_dir):
        if pdb_file.endswith(".pdb"):
            pdb_file_path = os.path.join(protein_dir, pdb_file)
            process_protein(pdb_file_path, redock_dir, mgltool_path, autodock_gpu_path, karmadock_path, ledock_path)

if __name__ == "__main__":
    main()
