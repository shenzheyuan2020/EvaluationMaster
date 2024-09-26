import requests
import csv
import pandas as pd
import sys
import os

def download_alphafold_pdb(uniprot_id, save_dir):
    """Download PDB file from AlphaFold database given a UniProt ID."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        os.makedirs(save_dir, exist_ok=True)
        file_name = f"AF-{uniprot_id}-F1-model_v3.pdb"
        file_path = os.path.join(save_dir, file_name)
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"Downloaded AlphaFold PDB: {file_path}")
    else:
        print(f"Failed to download AlphaFold PDB for UniProt ID {uniprot_id}. Status Code: {response.status_code}. URL: {url}")

def download_xray_pdb(pdb_id, save_dir):
    """Download X-ray PDB file from RCSB database given a PDB ID."""
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)
    if response.status_code == 200:
        os.makedirs(save_dir, exist_ok=True)
        file_path = os.path.join(save_dir, f"{pdb_id}.pdb")
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"Downloaded X-ray PDB: {file_path}")
    else:
        print(f"Failed to download X-ray PDB for {pdb_id}")

def download_nmr_pdb(pdb_id, save_dir):
    """Download NMR PDB file from RCSB database given a PDB ID."""
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)
    if response.status_code == 200:
        os.makedirs(save_dir, exist_ok=True)
        file_path = os.path.join(save_dir, f"{pdb_id}.pdb")
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"Downloaded NMR PDB: {file_path}")
    else:
        print(f"Failed to download NMR PDB for {pdb_id}")

def main(uniprot_id, name, resolution_threshold, method_filter, save_path):
    pdb_data = []

    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.txt'
    response = requests.get(url)

    if response.status_code == 200:
        data = response.text
        
        data_file_name = f'{name}_data.txt'
        with open(data_file_name, 'w') as f:
            f.write(data)
        print(f"UniProt data Downloaded: {data_file_name}")

        for line in data.splitlines():
            if line.startswith("DR   PDB;"):
                parts = line.split(";")
                pdb_id = parts[1].strip()  
                method = parts[2].strip()  
                resolution = parts[3].strip()  
                pdb_data.append([pdb_id, resolution, method])  

        csv_file_path = f'{save_path}/{name}_pdb_info.csv'
        with open(csv_file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['PDB ID', 'Resolution', 'Method'])  
            writer.writerows(pdb_data)  
        
        print(f"pdb info saved: {csv_file_path}")
        
        # Create a DataFrame from the PDB data
        df = pd.DataFrame(pdb_data, columns=['PDB ID', 'Resolution', 'Method'])
        
        df['Resolution'] = df['Resolution'].str.replace(' A', '').replace('-', '0.0')  
        df['Resolution'] = pd.to_numeric(df['Resolution'], errors='coerce')
    
        for index, row in df.iterrows():
            pdb_id = row['PDB ID']
            method = row['Method'].strip()

            if method_filter.lower() == 'xray' and row['Resolution'] < resolution_threshold:
                download_xray_pdb(pdb_id, os.path.join(save_path, f"{name}_pdb"))
            elif method_filter.lower() == 'nmr' and method.lower() == 'nmr':  
                download_nmr_pdb(pdb_id, os.path.join(save_path, f"{name}_pdb"))
            # elif method_filter.lower() == 'alphafold':
            #     download_alphafold_pdb(uniprot_id, os.path.join(save_path, f"{name}_pdb"))
                break

        download_alphafold_pdb(uniprot_id, os.path.join(save_path, f"{name}_pdb"))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5])
