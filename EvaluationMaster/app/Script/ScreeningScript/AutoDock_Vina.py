import csv
import os
import shutil
import subprocess
import glob
import os
import sys
from vina import Vina
v = Vina(sf_name='vinardo')
import openpyxl
from openpyxl import Workbook



def run_command(command):
    subprocess.run(command, shell=True)

def prepare_receptor(base_dir , mol_name, tool_path):
    """Prepare a specific receptor"""
    os.chdir (base_dir)
    pdbname = f"{base_dir}/{mol_name}"
    prepare_receptor_cmd = f"{tool_path}/bin/pythonsh {tool_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r {pdbname}.pdb -o {pdbname}.pdbqt -A checkhydrogens -U nphs_lps_waters"
    run_command(prepare_receptor_cmd)


def run_dock_command(base_dir, mol_name, center, lig_path):
    pdbname = f"{base_dir}/{mol_name}"
    os.chdir(base_dir)
    v = Vina(sf_name='vinardo')
    v.set_receptor(f"{pdbname}.pdbqt")
    pdbqt_files = [f for f in os.listdir(lig_path) if f.endswith(".pdbqt")]

    results = []  # List to store name and score for sorting

    for pdbqt_file in pdbqt_files:
        lig_name = pdbqt_file.split(".")[0]
        lig_parm = os.path.join(base_dir, lig_name)
        file_path = os.path.join(lig_path, pdbqt_file)
        try:
            dock_command = f"vina --receptor {pdbname}.pdbqt --ligand {file_path} --center_x {center[0]} --center_y {center[1]} --center_z {center[2]} --size_x 30 --size_y 30 --size_z 30 --exhaustiveness=32 --out {lig_parm}.pdbqt"
            run_command(dock_command)

            with open(f"{lig_parm}.pdbqt", 'r') as output_file:
                lines = output_file.readlines()
                if len(lines) >= 2:
                    score_line = lines[1]
                    score = float(score_line.split()[3])
                    results.append([lig_name, score])
        except Exception as e:
            print(f"Error occurred while processing {pdbqt_file}: {e}")
            continue

    results.sort(key=lambda x: x[1])  # Sort by score in ascending order

    with open(f"{base_dir}/results.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Name", "Score"])
        for result in results:
            writer.writerow(result)

    # Add smiles column from ref_standard.csv
    add_smiles_column(os.path.join(base_dir, '..', '..', 'ref_standard.csv'), f"{base_dir}/results.csv", f"{base_dir}/../results.csv")


def add_smiles_column(ref_csv_path, docking_results_path, result_csv_path):
    # Build a lookup dictionary for smiles based on Name
    smiles_lookup = {}
    with open(ref_csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            smiles_lookup[row["Name"]] = row["smiles"]

    # Read docking results and add smiles column
    with open(docking_results_path, newline='') as csvfile, open(result_csv_path, 'w', newline='') as resultfile:
        reader = csv.reader(csvfile)
        writer = csv.writer(resultfile)
        header = next(reader)
        writer.writerow(["smiles"] + header)  # Add smiles column to header
        for row in reader:
            name = row[0]
            smiles = smiles_lookup.get(name, "")
            writer.writerow([smiles] + row)



def process_molecule(mol_name, center, base_dir, lig_path, tool_path):
    print(f"Processing {mol_name} with Geometric Center: {center}")
    prepare_receptor(base_dir, mol_name, tool_path)
    #Generate_Grid(base_dir, mol_name, center,lig_path)
    run_dock_command(base_dir, mol_name,  center, lig_path)
    # Note: The Glide running steps have been removed as per your request

def process_vina_csv(csv_file_path, lig_path, tool_path, save_path):
    (base_dir, base_name) = os.path.split(csv_file_path)
    #current_path = os.getcwd()
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skip the header row
        for row in reader:
            mol_name, x, y, z = row[0], float(row[1]), float(row[2]), float(row[3])
            mol_name_path = mol_name + "_Vina"
            os.makedirs(f"{save_path}/{mol_name_path}", exist_ok=True)
            shutil.copy(f"{base_dir}/{mol_name}.pdb",f"{save_path}/{mol_name_path}/{mol_name}.pdb")
            process_molecule(mol_name, (x, y, z), f"{save_path}/{mol_name_path}", lig_path, tool_path )

def main():

    csv_file= sys.argv[1]
    lig_path =sys.argv[2]
    tool_path = sys.argv[3]
    save_path = sys.argv[4]

    process_vina_csv(csv_file, lig_path, tool_path, save_path)


if __name__ == "__main__":
    main()

