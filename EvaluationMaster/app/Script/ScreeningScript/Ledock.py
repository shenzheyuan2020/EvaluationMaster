
import csv
import os
import glob
import threading
import subprocess
import sys
import shutil
from openpyxl import Workbook

def extract_score(dok_file):
    with open(dok_file, 'r') as file:
        lines = file.readlines()
        if len(lines) > 1:
            parts = lines[1].split()
            score = parts[-2]  # 倒数第二个元素是得分
            return score
    return None

def write_scores_to_csv(dok_files, directory):
    csv_path = os.path.join(directory, "results.csv")
    scores = []
    for dok_file in dok_files:
        name = os.path.splitext(os.path.basename(dok_file))[0]  # Remove the file extension
        score = extract_score(dok_file)
        if score:
            scores.append((name, float(score)))  # Convert score to float for sorting

    # Sort the scores in ascending order
    scores.sort(key=lambda x: x[1])

    # Write the sorted results to CSV
    with open(csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Name', 'Score'])
        for name, score in scores:
            writer.writerow([name, score])

    print(f"CSV file with docking results saved to {csv_path}")
    add_smiles_and_save_results(directory)

def add_smiles_and_save_results(directory):
    # 定位ref_standard.csv文件
    ref_standard_path = os.path.join(directory, '..', '..', 'ref_standard.csv')
    results_path = os.path.join(directory, '..', 'results.csv')

    # 读取ref_standard.csv中的SMILES信息
    smiles_dict = {}
    with open(ref_standard_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            smiles_dict[row['Name']] = row['smiles']

    # 读取当前目录中的docking_results.csv文件
    docking_results_path = os.path.join(directory, 'results.csv')
    updated_rows = []
    with open(docking_results_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row['Name']
            smiles = smiles_dict.get(name, '')  # 获取SMILES，如果不存在则返回空字符串
            updated_rows.append([smiles, name, row['Score']])  # SMILES作为第一列

    # 将更新后的行写入新的CSV文件
    with open(results_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['smiles', 'Name', 'Score'])  # 写入头部
        writer.writerows(updated_rows)

    print(f"Updated results saved to {results_path}")


def create_in_file(directory, mol_name, center, lig_path, thread_id, total_threads, ledock_path, csv_file_path):
    csv_file_name = os.path.split(csv_file_path)[0]
    x1 = center[0] - 10
    x2 = center[0] + 10
    y1 = center[1] - 10
    y2 = center[1] + 10
    z1 = center[2] - 10
    z2 = center[2] + 10
    content = f"""\
Receptor
{csv_file_name}/{mol_name}.pdb

RMSD
1.0

Binding pocket
{x1:.3f} {x2:.3f}
{y1:.3f} {y2:.3f}
{z1:.3f} {z2:.3f}

Number of binding poses
30

Ligands list
{lig_path}/ligands_{thread_id}

END
"""
    os.makedirs(directory, exist_ok=True)
    dock_file = f"dock_{thread_id}.in"
    with open(os.path.join(directory, dock_file), "w") as f:
        f.write(content)

    # Splitting .mol2 files into batches
    ligands_list = sorted(glob.glob(os.path.join(lig_path, "*.mol2")))
    batch_size = len(ligands_list) // total_threads
    start_index = thread_id * batch_size
    end_index = None if thread_id == total_threads - 1 else start_index + batch_size
    batch_ligands = ligands_list[start_index:end_index]

    ligands_file = f"ligands_{thread_id}"
    with open(os.path.join(lig_path, ligands_file), "w") as file:
        for ligand in batch_ligands:
            file.write(os.path.abspath(ligand) + "\n")
    print(f"Ligands list {ligands_file} created in {lig_path}")

    # Running LEDock in a separate thread
    shell_command = f"{ledock_path}/ledock {os.path.join(directory, dock_file)}"
    try:
        os.chdir(directory)
        subprocess.run(shell_command, shell=True, check=True)
        print(f"LEDock for thread {thread_id} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"LEDock for thread {thread_id} failed with error code {e.returncode}.")

def process_molecule(mol_name, center, base_dir, lig_path, thread_id, total_threads, ledock_path , csv_file_path):
    print(f"Processing {mol_name} with Geometric Center: {center} for thread {thread_id}")
    create_in_file(base_dir, mol_name, center, lig_path, thread_id, total_threads, ledock_path, csv_file_path)

def move_and_overwrite(dok_file, mol_name_path):
    target_file = os.path.join(mol_name_path, os.path.basename(dok_file))
    if os.path.exists(target_file):
        os.remove(target_file)  # 如果目标文件已存在，先删除它
    shutil.move(dok_file, mol_name_path)  # 然后移动文件


def process_ledock_csv(csv_file_path, lig_path, total_threads, save_path, ledock_path):
    (base_dir, base_name) = os.path.split(csv_file_path)
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skip the header row
        for row in reader:
            mol_name, x, y, z = row[0], float(row[1]), float(row[2]), float(row[3])
            mol_name_path = f"{save_path}/{mol_name}_ledock"
            os.makedirs(mol_name_path, exist_ok=True)
            threads = []
            for thread_id in range(total_threads):
                thread = threading.Thread(target=process_molecule, args=(mol_name, (x, y, z), mol_name_path, lig_path, thread_id, total_threads, ledock_path, csv_file_path))
                threads.append(thread)
                thread.start()
            for thread in threads:
                thread.join()

            # Move .dok files to the corresponding mol_name directory after all threads have completed
            dok_files = glob.glob(os.path.join(lig_path, "*.dok"))
            for dok_file in dok_files:
                move_and_overwrite(dok_file, mol_name_path)
            print(f"All .dok files moved to {mol_name_path}")
            all_dok_files = glob.glob(os.path.join(mol_name_path, "*.dok"))
            write_scores_to_csv(all_dok_files, mol_name_path)
    # After all molecules have been processed and all threads have finished, collect all .dok files
    # This assumes dok files are located directly within the mol_name_path directories

def main():
    total_threads = int(sys.argv[1])
    csv_file_path = sys.argv[2]
    lig_path = sys.argv[3]
    save_path = sys.argv[4]
    ledock_path = sys.argv[5]
    process_ledock_csv(csv_file_path, lig_path, total_threads, save_path, ledock_path)

if __name__ == "__main__":
    main()
