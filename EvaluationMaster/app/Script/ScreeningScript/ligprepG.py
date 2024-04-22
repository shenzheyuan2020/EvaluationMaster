import os
import subprocess
import shutil
from file_conversion_script import convert_file_to_csv, process_csv_file  # 从第一个脚本导入转换函数
import sys

def generate_inp_file(csv_filename):
    """
    Generates the content of an INP file based on the given CSV filename.
    """
    inp_content = f"""INPUT_FILE_NAME   {csv_filename}
OUT_MAE   ligands.maegz
MAX_ATOMS   500
FORCE_FIELD   16
EPIK   yes
DETERMINE_CHIRALITIES   no
IGNORE_CHIRALITIES   no
NUM_STEREOISOMERS   32
"""
    return inp_content

def generate_sh_file(inp_filename, Schro_dir, njobs=12):
    """
    Generates the content of a SH file based on the given INP filename.
    """
    sh_content = f""" "{Schro_dir}/ligprep" -inp {inp_filename} -HOST localhost:{njobs} -NJOBS {njobs} -TMPLAUNCHDIR -WAIT
"""
    return sh_content

def process_molecules(njobs, input_file, save_dir, Schro_dir):
    # Check file extension and convert if necessary
    file_ext = os.path.splitext(input_file)[1].lower()
    if file_ext == '.smi':
        csv_filename = input_file
    elif file_ext == '.csv':
        csv_filename = input_file
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    # Ensure the working directory exists
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Directly use save_dir as the output directory
    csv_basename = os.path.basename(csv_filename)
    csv_destination = os.path.join(save_dir, csv_basename)
    shutil.copy(csv_filename, csv_destination)
    print(f"CSV file copied to: {csv_destination}")

    # Generate and save the INP file
    inp_filepath = os.path.join(save_dir, os.path.splitext(csv_basename)[0] + '.inp')
    inp_content = generate_inp_file(csv_basename)
    with open(inp_filepath, 'w') as inp_file:
        inp_file.write(inp_content)
    print(f"INP file created at: {inp_filepath}")

    # Generate and save the SH file
    sh_filepath = os.path.join(save_dir, os.path.splitext(csv_basename)[0] + '.sh')
    sh_content = generate_sh_file(os.path.basename(inp_filepath), Schro_dir, njobs)
    with open(sh_filepath, 'w') as sh_file:
        sh_file.write(sh_content)
    print(f"SH file created at: {sh_filepath}")

    # Change the working directory to save_dir
    os.chdir(save_dir)
    print(f"Changed working directory to: {os.getcwd()}")

    # Execute the SH file
    try:
        sh_file_name = os.path.basename(sh_filepath)
        subprocess.run(['bash', sh_file_name])
        print(f"Executed SH file: {sh_file_name}")
    except Exception as e:
        print(f"Error executing SH file: {e}")

# def main(input_file, working_directory):
#     current_path = os.getcwd()
#     process_molecules(input_file, f'{current_path}/{working_directory}')

if __name__ == "__main__":

    njobs = sys.argv[1]
    input_file = sys.argv[2]
    save_dir = sys.argv[3]
    Schro_dir = sys.argv[4]
    process_molecules(njobs, input_file, save_dir, Schro_dir)








    # # Example usage when running this script directly
    # input_file = 'working_place/DDR1_clus_decoys-selected.smi'  # 或者 .csv 文件
    # working_directory = 'working_place/'
    # main(input_file, working_directory)

# Note: Ensure file_conversion_script.py is in the same directory or in the PYTHONPATH
