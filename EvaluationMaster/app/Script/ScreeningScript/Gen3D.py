###ligprep4has absolute path


import os
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
def convert_sdf_to_mol2(sdf_path, mol2_path):
    """
    使用OpenBabel将SDF文件转换为MOL2文件。

    :param sdf_path: SDF文件的路径。
    :param mol2_path: 输出的MOL2文件的路径。
    """
    try:
        subprocess.run(["obabel", sdf_path, "-O", mol2_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error in converting {sdf_path} to MOL2: {e}")

def convert_sdf_to_pdbqt(mol2_path, pdbqt_path, mgltool_path):
    """
    使用MGLTOOLS将mol2文件转换为PDBQT文件。

    :param mol2_path: MOL2文件的路径。
    :param pdbqt_path: 输出的PDBQT文件的路径。
    :param mgltool_path: MGLTOOLS的路径。
    """
    try:
        subprocess.run([
                            f"{mgltool_path}/bin/pythonsh", 
                            f"{mgltool_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py", 
                            "-l", mol2_path, 
                            "-o", pdbqt_path
                        ], check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error in converting {mol2_path} to PDBQT: {e}")


def optimize_and_select_lowest_energy_conformer(mol):
    """
    对分子的所有构象进行优化，并选择能量最低的构象。

    :param mol: RDKit分子对象。
    :return: 能量最低的构象ID，如果无法优化则返回None。
    """
    try:
        AllChem.MMFFSanitizeMolecule(mol)
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol)

        lowest_energy = float('inf')
        lowest_energy_confId = None

        for confId in range(mol.GetNumConformers()):
            ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=confId)
            if ff is None:  # 检查力场是否有效
                continue
            energy = ff.CalcEnergy()
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_energy_confId = confId

        return lowest_energy_confId
    except Exception as e:
        print(f"Error during optimization: {e}")
        return None
    

def generate_3d_structures(csv_file, save_dir, mgltool_path):
    """
    从CSV文件中读取SMILES字符串和名称，生成它们的3D结构，并将它们保存在不同的文件夹中。

    :param csv_file: 包含SMILES字符串和名称的CSV文件的路径。
    """
    df = pd.read_csv(csv_file)

    # 获取CSV文件所在目录并创建子目录
    # directory = os.path.dirname(csv_file)
    # csv_basename = os.path.basename(csv_file)
    # subdir_name = os.path.splitext(csv_basename)[0]
    # file_dir= subdir_name.split("_")[0]
    # subdir_path = os.path.join(save_dir, file_dir)
    subdir_path = save_dir




    sdf_output_directory = os.path.join(subdir_path, 'SDF')
    mol2_output_directory = os.path.join(subdir_path, 'MOL2')
    pdbqt_output_directory = os.path.join(subdir_path, 'PDBQT')

    for d in [sdf_output_directory, mol2_output_directory, pdbqt_output_directory]:
        if not os.path.exists(d):
            os.makedirs(d)

    for index, row in df.iterrows():
        smi = row['smiles']
        name = row['Name']  # 直接使用CSV文件中的Name作为basename
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        m2 = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(m2, numConfs=10, numThreads=4)

        # 选择能量最低的构象
        lowest_energy_confId = optimize_and_select_lowest_energy_conformer(m2)
        if lowest_energy_confId is None:
            continue

        # 保存为SDF格式并转换为MOL2格式
        output_sdf = os.path.join(sdf_output_directory, f'{name}.sdf')
        output_mol2 = os.path.join(mol2_output_directory, f'{name}.mol2')
        output_pdbqt= os.path.join(pdbqt_output_directory, f'{name}.pdbqt')
        sdf_writer = Chem.SDWriter(output_sdf)
        sdf_writer.write(m2, confId=lowest_energy_confId)
        sdf_writer.close()

        # 转换SDF到MOL2
        convert_sdf_to_mol2(output_sdf, output_mol2)
        # 转换MOL2到pdbqt
        os.chdir(mol2_output_directory)
        convert_sdf_to_pdbqt(output_mol2, output_pdbqt, mgltool_path)

def main():
    # 替换以下路径为您的实际文件路径
    csv_file_path = sys.argv[1]
    save_dir = sys.argv[2]
    mgltool_path = sys.argv[3]
    # csv_file_path = 'working_place/lig/DDR1_clus_standard.csv'


    generate_3d_structures(csv_file_path, save_dir, mgltool_path)
if __name__ == "__main__":
    main()







