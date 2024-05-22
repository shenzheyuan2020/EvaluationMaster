# ###ligprep4has absolute path


# import os
# import subprocess
# import pandas as pd
# from rdkit import Chem
# from rdkit.Chem import AllChem
# import sys
# def convert_sdf_to_mol2(sdf_path, mol2_path):
#     """
#     使用OpenBabel将SDF文件转换为MOL2文件。

#     :param sdf_path: SDF文件的路径。
#     :param mol2_path: 输出的MOL2文件的路径。
#     """
#     try:
#         subprocess.run(["obabel", sdf_path, "-O", mol2_path], check=True)
#     except subprocess.CalledProcessError as e:
#         print(f"Error in converting {sdf_path} to MOL2: {e}")

# def convert_sdf_to_pdbqt(mol2_path, pdbqt_path, mgltool_path):
#     """
#     使用MGLTOOLS将mol2文件转换为PDBQT文件。

#     :param mol2_path: MOL2文件的路径。
#     :param pdbqt_path: 输出的PDBQT文件的路径。
#     :param mgltool_path: MGLTOOLS的路径。
#     """
#     try:
#         subprocess.run([
#                             f"{mgltool_path}/bin/pythonsh", 
#                             f"{mgltool_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py", 
#                             "-l", mol2_path, 
#                             "-o", pdbqt_path
#                         ], check=True)

#     except subprocess.CalledProcessError as e:
#         print(f"Error in converting {mol2_path} to PDBQT: {e}")


# def optimize_and_select_lowest_energy_conformer(mol):
#     """
#     对分子的所有构象进行优化，并选择能量最低的构象。

#     :param mol: RDKit分子对象。
#     :return: 能量最低的构象ID，如果无法优化则返回None。
#     """
#     try:
#         AllChem.MMFFSanitizeMolecule(mol)
#         mmff_props = AllChem.MMFFGetMoleculeProperties(mol)

#         lowest_energy = float('inf')
#         lowest_energy_confId = None

#         for confId in range(mol.GetNumConformers()):
#             ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=confId)
#             if ff is None:  # 检查力场是否有效
#                 continue
#             energy = ff.CalcEnergy()
#             if energy < lowest_energy:
#                 lowest_energy = energy
#                 lowest_energy_confId = confId

#         return lowest_energy_confId
#     except Exception as e:
#         print(f"Error during optimization: {e}")
#         return None
    

# def generate_3d_structures(csv_file, lig_start_name, save_dir, mgltool_path):
#     """
#     从CSV文件中读取SMILES字符串和名称，生成它们的3D结构，并将它们保存在不同的文件夹中。

#     :param csv_file: 包含SMILES字符串和名称的CSV文件的路径。
#     """
#     df = pd.read_csv(csv_file)

#     # 获取CSV文件所在目录并创建子目录
#     # directory = os.path.dirname(csv_file)
#     csv_basename = os.path.basename(csv_file)
#     subdir_name = os.path.splitext(csv_basename)[0]
#     file_dir= subdir_name.split("_")[0]
#     subdir_path = os.path.join(save_dir, file_dir)




#     sdf_output_directory = os.path.join(subdir_path, 'SDF')
#     mol2_output_directory = os.path.join(subdir_path, 'MOL2')
#     pdbqt_output_directory = os.path.join(subdir_path, 'PDBQT')

#     for d in [sdf_output_directory, mol2_output_directory, pdbqt_output_directory]:
#         if not os.path.exists(d):
#             os.makedirs(d)

#     for index, row in df.iterrows():
#         smi = row['smiles']
#         if lig_start_name == "original_name":
#             name = row['Name']
#         else:
#             name = f"{lig_start_name}_{index + 1}"
#         mol = Chem.MolFromSmiles(smi)
#         if mol is None:
#             continue

#         m2 = Chem.AddHs(mol)
#         AllChem.EmbedMultipleConfs(m2, numConfs=10, numThreads=4)

#         # 选择能量最低的构象
#         lowest_energy_confId = optimize_and_select_lowest_energy_conformer(m2)
#         if lowest_energy_confId is None:
#             continue

#         # 保存为SDF格式并转换为MOL2格式
#         output_sdf = os.path.join(sdf_output_directory, f'{name}.sdf')
#         output_mol2 = os.path.join(mol2_output_directory, f'{name}.mol2')
#         output_pdbqt= os.path.join(pdbqt_output_directory, f'{name}.pdbqt')
#         sdf_writer = Chem.SDWriter(output_sdf)
#         sdf_writer.write(m2, confId=lowest_energy_confId)
#         sdf_writer.close()

#         # 转换SDF到MOL2
#         convert_sdf_to_mol2(output_sdf, output_mol2)
#         # 转换MOL2到pdbqt
#         os.chdir(mol2_output_directory)
#         convert_sdf_to_pdbqt(output_mol2, output_pdbqt, mgltool_path)

# def main():
#     # 替换以下路径为您的实际文件路径
#     lig_start_name = sys.argv[1]
#     csv_file_path = sys.argv[2]
#     save_dir = sys.argv[3]
#     mgltool_path = sys.argv[4]
#     # csv_file_path = 'working_place/lig/DDR1_clus_standard.csv'


#     generate_3d_structures(csv_file_path, lig_start_name, save_dir, mgltool_path)
# if __name__ == "__main__":
#     main()





import os
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def convert_sdf_to_mol2(sdf_path, mol2_path):
    """
    使用OpenBabel将SDF文件转换为MOL2文件。

    :param sdf_path: SDF文件的路径。
    :param mol2_path: 输出的MOL2文件的路径。
    """
    try:
        subprocess.run(["obabel", sdf_path, "-O", mol2_path], check=True)
        logging.info(f"Successfully converted {sdf_path} to {mol2_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in converting {sdf_path} to MOL2: {e}")

def convert_sdf_to_pdbqt(mol2_path, pdbqt_path, mgltool_path):
    """
    使用MGLTOOLS将mol2文件转换为PDBQT文件。

    :param mol2_path: MOL2文件的路径。
    :param pdbqt_path: 输出的PDBQT文件的路径。
    :param mgltool_path: MGLTOOLS的路径。
    """
    try:
        subprocess.run([
            os.path.join(mgltool_path, 'bin', 'pythonsh'),
            os.path.join(mgltool_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_ligand4.py'),
            '-l', mol2_path,
            '-o', pdbqt_path
        ], check=True)
        logging.info(f"Successfully converted {mol2_path} to {pdbqt_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in converting {mol2_path} to PDBQT: {e}")

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
        logging.error(f"Error during optimization: {e}")
        return None

def process_molecule(index, row, lig_start_name, sdf_output_directory, mol2_output_directory, pdbqt_output_directory, mgltool_path):
    smi = row['smiles']
    if lig_start_name == "original_name":
        name = row['Name'] if pd.notnull(row['Name']) else f"unnamed_{index + 1}"
    else:
        name = f"{lig_start_name}_{index + 1}"
    
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logging.warning(f"Failed to generate molecule from SMILES: {smi}")
            return

        m2 = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(m2, numConfs=10, numThreads=1)

        # 选择能量最低的构象
        lowest_energy_confId = optimize_and_select_lowest_energy_conformer(m2)
        if lowest_energy_confId is None:
            logging.warning(f"Failed to optimize conformers for molecule: {name}")
            return

        # 保存为SDF格式并转换为MOL2格式
        output_sdf = os.path.join(sdf_output_directory, f'{name}.sdf')
        output_mol2 = os.path.join(mol2_output_directory, f'{name}.mol2')
        output_pdbqt = os.path.join(pdbqt_output_directory, f'{name}.pdbqt')
        
        sdf_writer = Chem.SDWriter(output_sdf)
        sdf_writer.write(m2, confId=lowest_energy_confId)
        sdf_writer.close()
        logging.info(f"Saved {name} as SDF")

        # 转换SDF到MOL2
        convert_sdf_to_mol2(output_sdf, output_mol2)
        # 转换MOL2到pdbqt
        os.chdir(mol2_output_directory)
        convert_sdf_to_pdbqt(output_mol2, output_pdbqt, mgltool_path)

    except Exception as e:
        logging.error(f"Error processing molecule {name}: {e}")

def generate_3d_structures(csv_file, lig_start_name, save_dir, mgltool_path, num_threads=1):
    """
    从CSV文件中读取SMILES字符串和名称，生成它们的3D结构，并将它们保存在不同的文件夹中。

    :param csv_file: 包含SMILES字符串和名称的CSV文件的路径。
    :param num_threads: 并行线程的数量，默认值为1。
    """
    df = pd.read_csv(csv_file)

    csv_basename = os.path.basename(csv_file)
    subdir_name = os.path.splitext(csv_basename)[0]
    file_dir = subdir_name.split("_")[0]
    subdir_path = os.path.join(save_dir, file_dir)

    sdf_output_directory = os.path.join(subdir_path, 'SDF')
    mol2_output_directory = os.path.join(subdir_path, 'MOL2')
    pdbqt_output_directory = os.path.join(subdir_path, 'PDBQT')

    for d in [sdf_output_directory, mol2_output_directory, pdbqt_output_directory]:
        os.makedirs(d, exist_ok=True)

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [
            executor.submit(
                process_molecule, index, row, lig_start_name, sdf_output_directory, mol2_output_directory, pdbqt_output_directory, mgltool_path
            ) for index, row in df.iterrows()
        ]
        
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Error in future: {e}")

def main():
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Usage: python script.py <lig_start_name> <csv_file_path> <save_dir> <mgltool_path> [num_threads]")
        sys.exit(1)

    lig_start_name = sys.argv[1]
    csv_file_path = sys.argv[2]
    save_dir = sys.argv[3]
    mgltool_path = sys.argv[4]
    num_threads = int(sys.argv[5]) if len(sys.argv) == 6 else 1

    generate_3d_structures(csv_file_path, lig_start_name, save_dir, mgltool_path, num_threads)

if __name__ == "__main__":
    main()
