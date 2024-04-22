import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import os
import sys

def bitvector_to_array(bitvector):
    bitlist = list(bitvector.ToBitString())
    return np.array([int(b) for b in bitlist], dtype=np.uint8)

def read_and_filter_data(file_path, pIC50threshold):
    data = pd.read_csv(file_path)
    return data[data['pIC50'] >= int(pIC50threshold)]

def calculate_fingerprints(filtered_data):
    filtered_data['Molecule'] = filtered_data['smiles'].apply(lambda x: Chem.MolFromSmiles(x))
    filtered_data['Fingerprint'] = filtered_data['Molecule'].apply(
        lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=512))
    return np.array([bitvector_to_array(fp) for fp in filtered_data['Fingerprint']])

def perform_clustering_and_save(filtered_data, fingerprints, file_path):
    distances = cdist(fingerprints, fingerprints, metric='jaccard')
    # cluster = AgglomerativeClustering(n_clusters=50, affinity='precomputed', linkage='average')
    cluster = AgglomerativeClustering(n_clusters=int(clus_num), affinity='deprecated', linkage='average')
    clusters = cluster.fit_predict(distances)
    filtered_data['Cluster'] = clusters
    top_molecules = filtered_data.loc[filtered_data.groupby('Cluster')['pIC50'].idxmax()]
    output_file_path = file_path.replace('.csv', '_clus.csv')
    top_molecules.to_csv(output_file_path)

def tsne_visualization_and_save(fingerprints, clusters, file_path):
    plt.rcParams.update({'font.size': 18})  # 可以根据需要调整这里的值
    tsne = TSNE(n_components=2, random_state=0)
    fingerprints_tsne = tsne.fit_transform(fingerprints)
    plt.figure(figsize=(10, 8))
    plt.scatter(fingerprints_tsne[:, 0], fingerprints_tsne[:, 1], c=clusters, cmap='viridis', s=50, alpha=0.3)
    plt.title('t-SNE of Molecular Fingerprints')
    plt.xlabel('t-SNE feature 1')
    plt.ylabel('t-SNE feature 2')
    plt.colorbar(label='Cluster ID')

    # 提高分辨率，例如设定为 300 DPI
    figpath = file_path.replace('.csv', '_tSNE.png')
    plt.savefig(figpath, dpi=1200)
    plt.close()
    return figpath

def process_excel_file(handle_file, pIC50threshold, clus_num):
    filtered_data = read_and_filter_data(handle_file, pIC50threshold)
    if len(filtered_data) < int(clus_num):
        print(f"Skipping file due to insufficient data points.")
        return
    else:
        fingerprints = calculate_fingerprints(filtered_data)
        perform_clustering_and_save(filtered_data, fingerprints, handle_file)
        figpath = tsne_visualization_and_save(fingerprints, filtered_data['Cluster'], handle_file)
        print(f"Processed and saved data and t-SNE plot for the file.")
    return figpath
if __name__ == "__main__":
    # 指定单个文件路径
    clus_num = sys.argv[1]
    pIC50threshold = sys.argv[2]
    handle_file  = sys.argv[3]


    process_excel_file(handle_file, pIC50threshold, clus_num)




    # file_path = 'usp7.csv'  # 替换为你的文件路径
    # process_excel_file(file_path)
