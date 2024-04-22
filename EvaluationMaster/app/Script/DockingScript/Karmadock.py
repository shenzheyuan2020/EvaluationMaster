import os
import pandas as pd
from tqdm import tqdm
import argparse
import sys
import numpy as np
import torch.nn.functional as F
import torch.nn as nn
import pandas as pd
import torch.optim
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess

def set_project_dir_and_imports(KarmaDock_dir):
    sys.path.append(KarmaDock_dir)
    from architecture.KarmaDock_architecture import KarmaDock
    from utils.post_processing import correct_pos
    from utils.pre_processing import get_pocket_pure
    from utils.fns import Early_stopper, set_random_seed as imported_set_random_seed
    from dataset.dataloader_obj import PassNoneDataLoader
    from dataset.graph_obj import get_mol2_xyz_from_cmd
    from dataset.graph_obj import VSTestGraphDataset_FlyReload_CSV
    # Returning all imported modules, functions, classes, including get_mol2_xyz_from_cmd
    return KarmaDock, correct_pos, get_pocket_pure, Early_stopper, imported_set_random_seed, PassNoneDataLoader, get_mol2_xyz_from_cmd, VSTestGraphDataset_FlyReload_CSV


def move_molecule_to_center_and_convert_to_mol2(out_dir, target_center):
    """
    生成一个分子，并将其移动到指定的中心坐标，然后保存为MOL2文件。
    分子的SMILES字符串和保存路径是固定的。
    
    参数:
    - target_center (tuple): 目标中心坐标，格式为(x, y, z)。
    """
    # 固定的SMILES字符串
    smiles = 'C=CC(=O)N1CCC[C@H](C1)N2C3=C(C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=NC=N3)N'
    
    # 获取当前工作目录并构造保存路径
    save_path = os.path.join(f'{out_dir}', "tmp", f"test")

    # 确保保存目录存在
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    
    # 从SMILES字符串创建分子，并添加氢原子
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    # 生成3D构象
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    
    # 获取分子的构象
    conf = mol.GetConformer()
    
    # 计算分子当前的几何中心
    center = AllChem.ComputeCentroid(conf)
    
    # 计算移动向量
    move_vector = (target_center[0] - center.x, target_center[1] - center.y, target_center[2] - center.z)
    
    # 根据移动向量调整所有原子的位置
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + move_vector[0], pos.y + move_vector[1], pos.z + move_vector[2]))
    
    # 保存为临时SDF文件以便转换
    temp_sdf_path = save_path + ".sdf"
    with Chem.SDWriter(temp_sdf_path) as writer:
        writer.write(mol)
    
    # 使用subprocess调用Open Babel进行格式转换
    temp_mol2_path = save_path + ".mol2"
    command = f'obabel -isdf {temp_sdf_path} -omol2 -O {temp_mol2_path}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # 检查命令是否成功执行并清理临时SDF文件
    if result.returncode == 0:
        print("Conversion successful.")
        os.remove(temp_sdf_path)
    else:
        print("Conversion failed.")
        print("Error:", result.stderr)
    return temp_mol2_path

def run_virtual_screening(protein_file, mol_name, crystal_ligand_file,  pocket_center, lig_csv, out_dir, score_threshold, KarmaDock_dir, VSTestGraphDataset_FlyReload_CSV):
    ligand_csv = f'{lig_csv}'
    graph_dir = os.path.join(f'{out_dir}/{mol_name}_graph/')
    os.makedirs(graph_dir, exist_ok=True)
    pocket_file = protein_file.replace('.pdb', '_pocket.pdb')
    batch_size = 64
    # get pocket center
    cry_lig_pos = get_mol2_xyz_from_cmd(crystal_ligand_file)
    # get pocket 
    if not os.path.exists(pocket_file):
        get_pocket_pure(protein_file, cry_lig_pos, pocket_file, size=12)
    # test 
    # Initialize test dataset, automatically generates graphs if needed
    test_dataset = VSTestGraphDataset_FlyReload_CSV(pocket_file, ligand_csv, pocket_center, graph_dir)
    test_dataset.generate_graphs(n_job=-1)  # Always ensure graphs are generated

    print('# virtual screening')
    vs_dir = os.path.join(f'{out_dir}/{mol_name}_KarmaDock/')
    os.makedirs(f'{out_dir}/{mol_name}_KarmaDock/', exist_ok=True)
    device_id = 0
    my_device = f'cuda:{device_id}' if torch.cuda.is_available() else 'cpu'
    
    model = KarmaDock()
    model = nn.DataParallel(model, device_ids=[device_id], output_device=device_id).to(my_device)
    model_file = f'{KarmaDock_dir}/trained_models/karmadock_screening.pkl'
    stopper = Early_stopper(model_file=model_file, mode='lower', patience=70)
    stopper.load_model(model_obj=model, my_device=my_device, strict=False)

    dst_csv = os.path.join(f'{out_dir}/{mol_name}_KarmaDock/', f'{mol_name}_score.csv')
    test_dataloader = PassNoneDataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=False, num_workers=0, follow_batch=[],  pin_memory=True)
    print(f'dataset num: {len(test_dataset)}')
    if len(test_dataset) > 0:
            model.eval()
            pki_scores_pred = torch.as_tensor([]).to(my_device)
            pki_scores_pred_ff = torch.as_tensor([]).to(my_device)
            pki_scores_pred_aligned = torch.as_tensor([]).to(my_device)
            pdb_ids = []
            with torch.no_grad():
                for idx, data in enumerate(tqdm(test_dataloader, desc='prediction')):
                    # try:
                    # to device
                    data = data.to(my_device)
                    batch_size = data['ligand'].batch[-1] + 1
                    # forward
                    pro_node_s, lig_node_s = model.module.encoding(data)
                    lig_pos, _, _ = model.module.docking(pro_node_s, lig_node_s, data, recycle_num=3)
                    mdn_score_pred = model.module.scoring(lig_s=lig_node_s, lig_pos=lig_pos, pro_s=pro_node_s, data=data,
                                                                                dist_threhold=5., batch_size=batch_size)
                    pki_scores_pred = torch.cat([pki_scores_pred, mdn_score_pred], dim=0)
                    pdb_ids.extend(data.pdb_id)
                    # # post processing
                    data.pos_preds = lig_pos
                    poses, _, _ = correct_pos(data, mask=mdn_score_pred <= score_threshold, out_dir=vs_dir, out_init=True, out_uncoorected=True, out_corrected=True)
                    ff_corrected_pos = torch.from_numpy(np.concatenate([i[0] for i in poses], axis=0)).to(my_device)
                    align_corrected_pos = torch.from_numpy(np.concatenate([i[1] for i in poses], axis=0)).to(my_device)
                    mdn_score_pred_ff_corrected = model.module.scoring(lig_s=lig_node_s, lig_pos=ff_corrected_pos, pro_s=pro_node_s, data=data,
                                                                            dist_threhold=5., batch_size=batch_size)
                    mdn_score_pred_align_corrected = model.module.scoring(lig_s=lig_node_s, lig_pos=align_corrected_pos, pro_s=pro_node_s, data=data,
                                                                        dist_threhold=5., batch_size=batch_size)
                    pki_scores_pred_ff = torch.cat([pki_scores_pred_ff, mdn_score_pred_ff_corrected], dim=0)
                    pki_scores_pred_aligned = torch.cat([pki_scores_pred_aligned, mdn_score_pred_align_corrected], dim=0)
                    # except:
                    #     continue
            pki_scores_pred = pki_scores_pred.view(-1).cpu().numpy().tolist()
            pki_scores_pred_ff = pki_scores_pred_ff.view(-1).cpu().numpy().tolist()
            pki_scores_pred_aligned = pki_scores_pred_aligned.view(-1).cpu().numpy().tolist()
            data = zip(pdb_ids, pki_scores_pred, pki_scores_pred_ff, pki_scores_pred_aligned) # pki_scores_pred_ff, pki_scores_pred_aligned
            columnds = ['Name', 'karma_score', 'karma_score_ff', 'karma_score_aligned']
            df = pd.DataFrame(data, columns=columnds)
            df.to_csv(dst_csv, index=False)


def main():
    KarmaDock_dir = sys.argv[1]
    multi_protein_csv = sys.argv[2]
    multi_protein_dir = sys.argv[3]
    lig_csv = sys.argv[4]
    out_dir = sys.argv[5]
    score_threshold = float(sys.argv[6])
    global KarmaDock, correct_pos, get_pocket_pure, Early_stopper, set_random_seed, PassNoneDataLoader, get_mol2_xyz_from_cmd, VSTestGraphDataset_FlyReload_CSV
    KarmaDock, correct_pos, get_pocket_pure, Early_stopper, set_random_seed, PassNoneDataLoader, get_mol2_xyz_from_cmd, VSTestGraphDataset_FlyReload_CSV = set_project_dir_and_imports(KarmaDock_dir)

    set_random_seed(2023)  
    proteins_dir = multi_protein_dir 
    data = pd.read_csv(multi_protein_csv)
    for index, row in data.iterrows():
        protein_file = os.path.join(proteins_dir, row['mol name'] + '.pdb')
        pocket_center = np.array([row['x'], row['y'], row['z']], dtype=float)
        mol_name = row['mol name']  
        print(f"Processing {mol_name}")
        target_center = pocket_center
        crystal_ligand_file = move_molecule_to_center_and_convert_to_mol2(out_dir, tuple(target_center))
        run_virtual_screening(protein_file, mol_name, crystal_ligand_file,  pocket_center, lig_csv, out_dir, score_threshold, KarmaDock_dir, VSTestGraphDataset_FlyReload_CSV)
        print(f"RESULT,{mol_name},{out_dir}")


if __name__ == '__main__':
    main()
