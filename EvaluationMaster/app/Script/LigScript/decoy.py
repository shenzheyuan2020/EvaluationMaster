import sys
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
# from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Chem import MolStandardize
import os
from collections import defaultdict
import numpy as np

from itertools import product
from joblib import Parallel, delayed
import re
from collections import defaultdict

# from IPython.display import clear_output
# IPythonConsole.ipython_useSVG = True


import pandas as pd
# Whether to use GPU for generating molecules with DeLinker
use_gpu = False
os.environ['CUDA_VISIBLE_DEVICES'] = '-1' if not use_gpu else '0'



def create_smi_from_csv(csv_file, smiles_column, output_smi_file):
    try:
        # 读取 CSV 文件
        df = pd.read_csv(csv_file)

        # 提取 SMILES 列
        smiles_data = df[smiles_column]

        # 将 SMILES 保存到 .smi 文件
        with open(output_smi_file, 'w') as file:
            for smile in smiles_data:
                file.write(smile + '\n')

        return output_smi_file
    except Exception as e:
        print(f"Error processing CSV file: {e}")
        return None

def smi_to_csv(smi_file, csv_file):
    try:
        # 读取 SMI 文件
        with open(smi_file, 'r') as file:
            lines = file.readlines()

        # 分析每行数据并去除第一列
        data = []
        decoy_counter = 1  # 用于生成递增的 decoy 名称
        for line in lines:
            parts = line.strip().split(' ')  # 使用空格分隔
            if len(parts) > 1:
                smiles = parts[1]  # 假设 SMILES 数据是每行的第二部分
                name = f"decoy_{decoy_counter}"  # 生成新的名称
                data.append([smiles, name])
                decoy_counter += 1
            else:
                # 如果没有两个部分，则跳过该行
                continue


        df = pd.DataFrame(data, columns=['smiles', 'Name'])


        df.to_csv(csv_file, index=False)

    except Exception as e:
        print(f"Error processing SMI file: {e}")



def process_smi_file(smi_file_path, generate_number, valid_number, tool_path):
    sys.path.append(tool_path)
    sys.path.append(f"{tool_path}/evaluation/")
    from DeepCoy import DenseGGNNChemModel
    from data.prepare_data import read_file, preprocess
    from select_and_evaluate_decoys import select_and_evaluate_decoys
    try:
        smi_file = os.path.basename(smi_file_path).split('.')[0]
        # Read and preprocess data
        raw_data = read_file(smi_file_path)
        preprocess(raw_data, "zinc", smi_file)

        # Set DeepCoy model arguments dynamically based on smi_file
        args = defaultdict(None)
        args['--dataset'] = 'zinc'
        args['--freeze-graph-model'] = False
        args['--restore'] =  f'{tool_path}/models/DeepCoy_DUDE_model_e09.pickle'

        args['--config'] = json.dumps({
            "generation": True,
            "batch_size": 1,
            "number_of_generation_per_valid": generate_number,
            "train_file": f"molecules_{smi_file}.json",
            "valid_file": f"molecules_{smi_file}.json",
            "output_name": f"{smi_file}_decoys.smi",
            "use_subgraph_freqs": False
        })

        # Initialize and train model
        model = DenseGGNNChemModel(args)
        model.train()
        del model  # Free memory

        # Evaluate decoys
        chosen_properties = "ALL"
        num_decoys_per_active = valid_number
        results = select_and_evaluate_decoys(smi_file + '_decoys.smi', 
                                            file_loc='./', output_loc='./', 
                                            dataset=chosen_properties, 
                                            num_cand_dec_per_act=num_decoys_per_active, 
                                            num_dec_per_act=num_decoys_per_active)
        smi_to_csv(smi_file + '_decoys-selected.smi', smi_file + '_decoys.csv')
        print(f"Processing {smi_file}")
        print("DOE score: \t\t\t%.3f" % results[8])
        print("Average Doppelganger score: \t%.3f" % results[10])
        print("Max Doppelganger score: \t%.3f" % results[11])
        print(f"成功处理 {smi_file}")
    except Exception as e:
        print(f"处理文件 {smi_file_path} 时发生错误: {e}")





def process_specific_file(csv_file_path, generate_number, valid_number, tool_path, save_dir):
    smiles_column='smiles'
    os.chdir(save_dir)
    try:

        base_name = os.path.splitext(os.path.basename(csv_file_path))[0]
        smi_file_name = f"{base_name}.smi"


        output_smi_file = os.path.join(os.path.dirname(csv_file_path), smi_file_name)
        

        smi_file_path = create_smi_from_csv(csv_file_path, smiles_column, output_smi_file)

        if smi_file_path:

            process_smi_file(smi_file_path, generate_number,valid_number, tool_path)
            print(f"成功处理 {smi_file_path}")
    except Exception as e:
        print(f"处理文件 {csv_file_path} 时发生错误: {e}")


def main():

    generate_number = int(sys.argv[1])
    valid_number = int(sys.argv[2])
    csv_file = sys.argv[3]
    tool_path  = sys.argv[4]
    save_dir = sys.argv[5]
    process_specific_file(csv_file, generate_number, valid_number, tool_path, save_dir)

if __name__ == '__main__':
    main()
