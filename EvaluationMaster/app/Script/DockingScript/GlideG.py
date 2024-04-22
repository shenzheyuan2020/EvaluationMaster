import csv
import subprocess
import shutil
import os
import sys
import pandas as pd

def create_in_file(directory, mol_name, center):
    content = f"""\
FORCEFIELD   OPLS4
GRID_CENTER   {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}
GRIDFILE   {mol_name}.zip
INNERBOX   10, 10, 10
OUTERBOX   30, 30, 30
RECEP_FILE   {mol_name}.maegz
"""
    os.chdir(directory)
    with open(f"{mol_name}.in", "w") as f:
        f.write(content)

def run_structconvert(csv_file_path, mol_name, base_dir, schrodinger_path, mmshare_path):
    schrodinger_path = os.environ.get('SCHRODINGER', f'{schrodinger_path}')  # Use environment variable if available
    pdb_path = os.path.split(csv_file_path)[0]
    pdb_path = f"{pdb_path}/{mol_name}"
    structconvert_command = [
        'run',
        f'{mmshare_path}/bin/Linux-x86_64/structconvert.py',
        f'{pdb_path}.pdb',
        f'{base_dir}/{mol_name}.maegz'
    ]
    os.chdir(base_dir)
    try:
        subprocess.run(structconvert_command, check=True)
        print("Structconvert completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Structconvert failed with error code {e.returncode}.")
        
def run_glide(mol_name, base_dir):
    glide_command = [
        'glide',
        f"{mol_name}.in",
        '-OVERWRITE',
        '-HOST', 'localhost',
        '-TMPLAUNCHDIR',
        '-WAIT'
    ]
    try:
        os.chdir(base_dir)
        subprocess.run(glide_command, check=True)
        print(f"Glide completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Glide failed with error code {e.returncode}.")

def create_new_in_file(directory, mol_name, lig_path, Mode):
    content = f"""\
FORCEFIELD   OPLS4
GRIDFILE   {mol_name}.zip
LIGANDFILE  {lig_path}
POSES_PER_LIG   30
POSTDOCK_NPOSE   30
PRECISION   {Mode}
"""
    os.chdir(directory)
    with open(f"{mol_name}_new.in", "w") as f:
        f.write(content)
        
def run_glide_adjust(mol_name, base_dir, threads_num):
    glide_command = [
        'glide',
        f"{mol_name}_new.in",
        '-OVERWRITE',
        '-adjust',
        '-HOST', f'localhost:{threads_num}',
        '-TMPLAUNCHDIR',
        '-WAIT'
    ]
    os.chdir(base_dir)
    try:
        subprocess.run(glide_command, check=True)
        print("Glide (adjust) completed successfully.")
        process_csv_after_glide_adjust(mol_name, base_dir)  # Process CSV after glide adjustment
    except subprocess.CalledProcessError as e:
        print(f"Glide (adjust) failed with error code {e.returncode}.")

def process_csv_after_glide_adjust(mol_name, base_dir):
    # 构建CSV文件路径
    csv_file_path = os.path.join(base_dir, f"{mol_name}_new.csv")
    # 读取CSV文件数据
    data = pd.read_csv(csv_file_path)
    # 只保留需要的两列，并重命名这两列为Name和Score
    data = data[['title', 'r_i_docking_score']].rename(columns={'title': 'Name', 'r_i_docking_score': 'Score'})
    # 先对数据进行排序：先按照Name进行排序，然后按照Score进行升序排序
    sorted_data = data.sort_values(by=['Name', 'Score'], ascending=[True, True])
    # 数据去重：在排序后的基础上，去除Name的重复项，仅保留每个Name对应的最小Score记录
    deduplicated_data = sorted_data.drop_duplicates('Name', keep='first')
    # 最后，按照Score列的大小进行最终排序
    final_sorted_data = deduplicated_data.sort_values(by='Score', ascending=True)
    # 构建Excel文件路径
    excel_file_path = os.path.join(base_dir, f"Docking_Results.xlsx")
    # 保存去重并排序后的数据到Excel
    final_sorted_data.to_excel(excel_file_path, index=False)
    # 打印保存信息
    print(f"Deduplicated and sorted data saved to {excel_file_path}")


# def move_and_convert_file(mol_name, base_dir):
#     src_file = f"{base_dir}/{mol_name}_new_pv.maegz"
#     dest_dir = f"{base_dir}/handle_result/{mol_name}/"
    
#     os.makedirs(dest_dir, exist_ok=True)
    
#     shutil.move(src_file, dest_dir)
#     print(f"Moved {src_file} to {dest_dir}")

#     structconvert_command = [
#         'run',
#         '/opt/schrodinger/2021-2/mmshare-v5.4/bin/Linux-x86_64/structconvert.py',
#         f"{dest_dir}/{mol_name}_new_pv.maegz",
#         f"{dest_dir}/{mol_name}.sdf"
#     ]
#     try:
#         subprocess.run(structconvert_command, check=True)
#         print("structconvert completed successfully.")
#     except subprocess.CalledProcessError as e:
#         print(f"structconvert failed with error code {e.returncode}.")

# def remove_duplicates_from_sdf(input_sdf_path, output_sdf_path):
#     seen_molecules = set()
#     buffer = []
#     with open(input_sdf_path, 'r') as f_in, open(output_sdf_path, 'w') as f_out:
#         for line in f_in:
#             if line.startswith("$$$$"):
#                 mol_name = buffer[0].strip()
#                 if mol_name not in seen_molecules:
#                     seen_molecules.add(mol_name)
#                     f_out.writelines(buffer)
#                     f_out.write("$$$$\n")
#                 buffer = []
#             else:
#                 buffer.append(line)

def process_molecule(csv_file_path, mol_name, center, base_dir, lig_path, schro_dir, mm_sahre_dir, Glide_mode, threads_num):
    print(f"Processing {mol_name} with Geometric Center: {center}")
    create_in_file(base_dir, mol_name, center)
    run_structconvert(csv_file_path, mol_name, base_dir, schro_dir, mm_sahre_dir)
    run_glide(mol_name, base_dir)
    create_new_in_file(base_dir, mol_name, lig_path, Glide_mode)
    run_glide_adjust(mol_name, base_dir, threads_num)
    # move_and_convert_file(mol_name, base_dir)

    # input_sdf_path = f"{base_dir}/handle_result/{mol_name}/{mol_name}.sdf"
    # output_sdf_path = f"{base_dir}/handle_result/{mol_name}/{mol_name}_no_duplicates.sdf"
    # remove_duplicates_from_sdf(input_sdf_path, output_sdf_path)

def process_glide_csv(csv_file_path, lig_path, schro_dir, mm_sahre_dir, save_dir, Glide_mode, threads_num):
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skip the header row
        for row in reader:
            mol_name, x, y, z = row[0], float(row[1]), float(row[2]), float(row[3])
            os.makedirs(os.path.join(save_dir, mol_name), exist_ok=True)
            process_molecule(csv_file_path, mol_name, (x, y, z), os.path.join(save_dir, mol_name), lig_path, schro_dir, mm_sahre_dir, Glide_mode, threads_num)

def main():

    # csv_file_path = "working_place/PDBhan/DDR1.csv"
    # lig_path = "working_place/DDR1/DDR1_clus_standard-out.maegz"
    csv_file_path = sys.argv[1]
    lig_file = sys.argv[2]
    schro_dir = sys.argv[3]
    mm_sahre_dir = sys.argv[4]
    save_dir = sys.argv[5]
    Glide_mode = sys.argv[6]
    threads_num = sys.argv[7]
    # lig_path = "working_place/DDR1/DDR1_clus_standard-out.maegz"
    process_glide_csv(csv_file_path, lig_file, schro_dir,  mm_sahre_dir, save_dir, Glide_mode, threads_num)



if __name__ == "__main__":
    main()