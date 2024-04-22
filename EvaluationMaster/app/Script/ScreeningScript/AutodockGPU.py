import csv
import os
import shutil
import subprocess
import glob
import os
import sys

def run_command(command):
    subprocess.run(command, shell=True)

def prepare_receptor(base_dir , mol_name,  mgltool_path):
    """Prepare a specific receptor"""
    pdbname = f"{base_dir}/{mol_name}"
    prepare_receptor_cmd = f"{mgltool_path}/bin/pythonsh {mgltool_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r {pdbname}.pdb -o {pdbname}.pdbqt -A checkhydrogens -U nphs_lps_waters"
    run_command(prepare_receptor_cmd)


def Generate_Grid(base_dir, mol_name, center, mgltool_path):
    os.chdir(base_dir)
    range = '50,50,50'
    pdbname = f"{base_dir}/{mol_name}"
    print(pdbname)
    center =f"{center[0]},{center[1]},{center[2]}"
    # create_gpf_file(pdbname, mol_name, center)
    prepare_gpf_cmd = f"{mgltool_path}/bin/pythonsh  {mgltool_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -r {pdbname}.pdbqt -o {pdbname}.gpf -p npts='{range}' -p gridcenter='{center}' -p ligand_types='HD,HS,C,A,N,NA,OA,F,SA,S,P,Cl,Br,I'"
    print(prepare_gpf_cmd)
    run_command(prepare_gpf_cmd)
    grid_generate_cmd = f"autogrid4 -p {pdbname}.gpf"
    run_command(grid_generate_cmd)

import subprocess
import os
import openpyxl
from openpyxl import Workbook

def extract_first_binding_energy(file_path):
    start_reading = False
    first_energy_score = None
    with open(file_path, 'r') as file:
        for line in file:
            if "_____|______|______|___________|_________|_________________|___________" in line:
                start_reading = True
                continue
            if start_reading and line.strip() and not first_energy_score:
                parts = line.split()
                if len(parts) > 3:
                    try:
                        first_energy_score = float(parts[3])
                        break
                    except ValueError:
                        continue
    return first_energy_score

# def run_dock_command(base_dir, mol_name, lig_path, Autdock_GPU_file, GPU_num, n_run):
#     pdbname = f"{base_dir}/{mol_name}"
#     os.chdir(base_dir)
#     pdbqt_files = [f for f in os.listdir(lig_path) if f.endswith(".pdbqt")]

#     # Initialize a workbook and sheet for Excel output
#     workbook = Workbook()
#     sheet = workbook.active
#     sheet.title = "Docking Results"
#     sheet.append(["Name", "Score"])  # Column titles

#     for pdbqt_file in pdbqt_files:
#         lig_name = pdbqt_file.split(".")[0]
#         file_path = os.path.join(lig_path, pdbqt_file)
#         lig_parm = os.path.join(base_dir, lig_name)
        
#         full_command = f"{Autdock_GPU_file} -ffile {pdbname}.maps.fld -lfile {file_path} -devnum {GPU_num} -resnam {lig_parm}  -nrun {n_run}"
#         try:
#             subprocess.run(full_command, shell=True, check=True)
#             print(f"处理文件 {pdbqt_file} 完成")
#             # Construct the path to the .dlg file
#             dlg_file_path = f"{lig_parm}.dlg"
#             # Extract the first binding energy
#             score = extract_first_binding_energy(dlg_file_path)
#             if score is not None:
#                 sheet.append([lig_name, score])
#         except subprocess.CalledProcessError as e:
#             print(f"处理文件 {pdbqt_file} 时出错: {e}")

#     # Save the workbook to an Excel file
#     excel_file_path = os.path.join(base_dir, "docking_results.xlsx")
#     workbook.save(filename=excel_file_path)
#     print(f"Excel file with docking results saved to {excel_file_path}")

def add_smiles_to_results(base_dir):
    results_file_name="results.csv" 
    ref_file_name="ref_standard.csv"
    # 构建路径到ref_standard.csv
    ref_file_path = os.path.join(base_dir, "..", "..", ref_file_name)
    results_file_path = os.path.join(base_dir, results_file_name)
    output_file_path = os.path.join(base_dir, "..", "results.csv")
    
    # 读取ref_standard.csv并构建一个以Name为键，smiles为值的字典
    smiles_dict = {}
    with open(ref_file_path, mode='r', newline='') as ref_file:
        reader = csv.DictReader(ref_file)
        for row in reader:
            smiles_dict[row["Name"]] = row["smiles"]
    
    # 读取results.csv，为每个Name添加smiles值
    results_with_smiles = []
    with open(results_file_path, mode='r', newline='') as results_file:
        reader = csv.DictReader(results_file)
        for row in reader:
            name = row["Name"]
            score = row["Score"]
            smiles = smiles_dict.get(name, "N/A")  # 如果找不到对应的smiles，则返回"N/A"
            results_with_smiles.append({"smiles": smiles, "Name": name, "Score": score})
    
    # 写入新的CSV文件，确保smiles列在第一列
    with open(output_file_path, mode='w', newline='') as output_file:
        fieldnames = ["smiles", "Name", "Score"]
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        for row in results_with_smiles:
            writer.writerow(row)
    
    print(f"CSV file with SMILES, docking results saved to {output_file_path}")

def run_dock_command(base_dir, mol_name, lig_path, Autdock_GPU_file, GPU_num, n_run):
    pdbname = f"{base_dir}/{mol_name}"
    os.chdir(base_dir)
    pdbqt_files = [f for f in os.listdir(lig_path) if f.endswith(".pdbqt")]

    results = []  # A list to store the docking results

    for pdbqt_file in pdbqt_files:
        lig_name = pdbqt_file.split(".")[0]
        file_path = os.path.join(lig_path, pdbqt_file)
        lig_parm = os.path.join(base_dir, lig_name)
        
        full_command = f"{Autdock_GPU_file} -ffile {pdbname}.maps.fld -lfile {file_path} -devnum {GPU_num} -resnam {lig_parm} -nrun {n_run}"
        try:
            subprocess.run(full_command, shell=True, check=True)
            print(f"处理文件 {pdbqt_file} 完成")
            dlg_file_path = f"{lig_parm}.dlg"
            score = extract_first_binding_energy(dlg_file_path)
            if score is not None:
                results.append((lig_name, score))
        except subprocess.CalledProcessError as e:
            print(f"处理文件 {pdbqt_file} 时出错: {e}")

    # Sort the results by score (the second item in the tuple)
    sorted_results = sorted(results, key=lambda x: x[1])

    # CSV file path
    csv_file_path = os.path.join(base_dir, "results.csv")
    
    # Writing to CSV
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Name", "Score"])  # Writing the headers
        for result in sorted_results:
            writer.writerow(result)  # Writing the data rows

    print(f"CSV file with docking results saved to {csv_file_path}")

    add_smiles_to_results(base_dir)

def process_molecule(mol_name, center, base_dir, lig_path,  mgltool_path, Autdock_GPU_file, GPU_num, n_run):
    print(f"Processing {mol_name} with Geometric Center: {center}")
    prepare_receptor(base_dir, mol_name, mgltool_path)
    Generate_Grid(base_dir, mol_name, center, mgltool_path)
    run_dock_command(base_dir, mol_name, lig_path,Autdock_GPU_file, GPU_num, n_run)
    # Note: The Glide running steps have been removed as per your request

def process_autodockgpu_csv(csv_file_path, lig_path, mgltool_path, Autdock_GPU_file, save_path, GPU_num, n_run):
    
    os.environ['CUDA_VISIBLE_DEVICES'] = str(GPU_num)
    (base_dir, base_name) = os.path.split(csv_file_path)
    #current_path = os.getcwd()
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skip the header row
        for row in reader:
            mol_name, x, y, z = row[0], float(row[1]), float(row[2]), float(row[3])
            mol_name_path = mol_name + "_Autodock_GPU"
            os.makedirs(f"{save_path}/{mol_name_path}", exist_ok=True)
            shutil.copy(f"{base_dir}/{mol_name}.pdb",f"{save_path}/{mol_name_path}/{mol_name}.pdb")
            process_molecule(mol_name, (x, y, z), f"{save_path}/{mol_name_path}", lig_path,  mgltool_path,  Autdock_GPU_file, GPU_num, n_run)

def main():
    # csv_file_path = "/home/hoo/Install/Deepevaluation/data_preparation/working_place/PDBhan/DDR1.csv"
    # lig_path = "/home/hoo/Install/Deepevaluation/data_preparation/working_place/lig/PDBQT/"
    # process_autodockgpu_csv(csv_file_path, lig_path)
    # #run_glide_adjust("working_place/PDBhan/", lig_path)


    csv_file= sys.argv[1]
    lig_path =sys.argv[2]
    mgltool_path = sys.argv[3]
    Autdock_GPU_file = sys.argv[4]
    save_path = sys.argv[5]
    GPU_num = int(sys.argv[6])-1
    n_run = sys.argv[7]
    process_autodockgpu_csv(csv_file, lig_path, mgltool_path, Autdock_GPU_file, save_path, GPU_num, n_run)


if __name__ == "__main__":
    main()
