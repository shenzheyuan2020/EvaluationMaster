#python predict_pocket.py /path/to/Support_software/DeepPocket/ /path/to/input/folder /path/to/coordinate.xlsx GPU_number
#example: python /home/hoo/Install/Evaluation_X/EvaluationMaster/app/Script/ProScript/predict_pocket.py /home/hoo/Install/Evaluation_X/Support_software/DeepPocket/ /home/hoo/Install/Evaluation_X/pdb_pocket_predict_test/ /home/hoo/Install/Evaluation_X/pdb_pocket_predict_test/coordinate.xlsx 0

from Bio.PDB import PDBParser, PDBIO, Select
import Bio
import pandas as pd
import os
import sys
import re
import torch.nn as nn
import torch
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score
import numpy as np
import importlib
import molgrid
import time
from skimage.morphology import binary_dilation
from skimage.morphology import cube
from skimage.morphology import closing
from skimage.segmentation import clear_border
from skimage.measure import label
import struct

#Load Deeppocket function
deep_pocket_path = sys.argv[1]
sys.path.append(deep_pocket_path)
from clean_pdb import clean_pdb
from get_centers import get_centers
from types_and_gninatyper import gninatype,create_types
from model import Model
from rank_pockets import test_model
from unet import Unet
from segment_pockets import test
import gc


def set_gpu(gpu_num):
    if torch.cuda.is_available():
        try:
            gpu_num = int(gpu_num)
            torch.cuda.set_device(gpu_num)
            print(f"Using GPU: {gpu_num}")
        except Exception as e:
            print(f"Specified GPU number is not available: {e}")
            sys.exit(1)
    else:
        print("No GPU available, will use CPU.")


def get_model_gmaker_eprovider(test_types, batch_size, model, checkpoint, dims=None):
    model.cuda()
    model.load_state_dict(checkpoint['model_state_dict'])
    eptest_large = molgrid.ExampleProvider(shuffle=False, stratify_receptor=False, labelpos=0, balanced=False, iteration_scheme=molgrid.IterationScheme.LargeEpoch, default_batch_size=batch_size)
    eptest_large.populate(test_types)
    if dims is None:
        gmaker = molgrid.GridMaker()
    else:
        gmaker = molgrid.GridMaker(dimension=dims)
    return model, gmaker, eptest_large
                
                
def process_pdb_files(input_folder, class_checkpoint_path, seg_checkpoint_path, rank=1, upsample='None', num_classes=1, threshold=0.5, mask_dist=3.5):
    pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')] 

    for pdb_file in pdb_files:
        try:
            protein_file = os.path.join(input_folder, pdb_file)
            protein_nowat_file = protein_file.replace('.pdb', '_nowat.pdb')

            # Clean PDB file
            clean_pdb(protein_file, protein_nowat_file)

            # Run fpocket
            os.system(f'fpocket -f {protein_nowat_file}')
            fpocket_dir = os.path.join(protein_nowat_file.replace('.pdb', '_out'), 'pockets')
            get_centers(fpocket_dir)
            barycenter_file = os.path.join(fpocket_dir, 'bary_centers.txt')

            # Process types and gninatyper
            protein_gninatype = gninatype(protein_nowat_file)
            class_types = create_types(barycenter_file, protein_gninatype)

            # Define seg_types and probs_types_file
            seg_types = class_types.replace('.types', '_ranked.types')
            probs_types_file = class_types.replace('.types', '_confidence.txt')

            # Rank pockets
            with torch.no_grad():
                class_model = Model()
                # class_checkpoint = torch.load(class_checkpoint_path)
                class_checkpoint = torch.load(class_checkpoint_path, map_location=torch.device('cpu'))

                with open(class_types, 'r') as file:
                    types_lines = file.readlines()
                batch_size = min(len(types_lines), 50)
                class_model, class_gmaker, class_eptest = get_model_gmaker_eprovider(class_types, batch_size, class_model, class_checkpoint)
                
                # Get ranking and probabilities
                class_labels, class_probs = test_model(class_model, class_eptest, class_gmaker, batch_size)
                sorted_zipped_lists = sorted(zip(class_probs[:len(types_lines)], types_lines), reverse=True)
                ranked_types = [element for _, element in sorted_zipped_lists]
                confidence_types = [str(element) for element, _ in sorted_zipped_lists]

            # Save rankings and probabilities to files
            with open(seg_types, 'w') as fout:
                fout.writelines(ranked_types)
            with open(probs_types_file, 'w') as fout:
                fout.writelines(f'{confidence}\n' for confidence in confidence_types)

        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")

        finally:
            # Clean up resources
            del class_model, class_checkpoint 
            gc.collect()
            torch.cuda.empty_cache()


def extract_coordinates(input_folder):
    data = []
    for item in os.listdir(input_folder):
        if os.path.isdir(os.path.join(input_folder, item)) and item.endswith('_nowat_out'):
            pockets_path = os.path.join(input_folder, item, 'pockets')
            bary_centers_file = os.path.join(pockets_path, 'bary_centers.txt')
            molecule_name = item.replace('_nowat_out', '')
            if os.path.isfile(bary_centers_file):
                with open(bary_centers_file, 'r') as file:
                    for line in file:
                        parts = line.strip().split()
                        if len(parts) == 4:
                            pocket_id, x, y, z = parts
                            data.append([f"{molecule_name}-{pocket_id}", x, y, z])
    return data


def save_coordinates_to_excel(data, excel_file_path):
    if not excel_file_path.endswith('.xlsx'):
        raise ValueError("Excel file path must end with '.xlsx'")

    pd.DataFrame(data, columns=['molecule_pocket', 'x', 'y', 'z']).to_excel(excel_file_path, index=False)


def main():
    if len(sys.argv) < 5:
        print("Sufficient path parameters are required, including the GPU number.")
        sys.exit(1)

    deep_pocket_path = sys.argv[1] 
    input_folder = sys.argv[2]
    excel_file_path = sys.argv[3]
    gpu_num = sys.argv[4]  

    # Set GPU
    set_gpu(gpu_num)

    # Define checkpoint paths
    class_checkpoint_path = os.path.join(deep_pocket_path, 'first_model_fold1_best_test_auc_85001.pth.tar')
    seg_checkpoint_path = os.path.join(deep_pocket_path, 'seg0_best_test_IOU_91.pth.tar')  
    
    # Process PDB files
    process_pdb_files(input_folder, class_checkpoint_path, seg_checkpoint_path) 

    # Extract and save coordinates
    coordinates = extract_coordinates(input_folder)
    save_coordinates_to_excel(coordinates, excel_file_path)


if __name__ == "__main__":
    main()
