import os
import sys
import requests

# 创建或确认 working_place/ 子目录的存在
# cwd = os.getcwd()
# working_dir = os.path.join(cwd, 'working_place')
# if not os.path.exists(working_dir):
#     os.makedirs(working_dir)

def download_pdb_file(pdb_id, name, save_dir):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)
    local_path = os.path.join(save_dir, name)
    if not os.path.exists(local_path):
        os.makedirs(local_path)

    file_path = os.path.join(save_dir, name, f'{pdb_id}.pdb')

    with open(file_path, 'w') as f:
        f.write(response.text)
    print(f"The PDB file has been saved to: {file_path}.")
    return file_path  # 返回下载文件的路径

if __name__ == "__main__":
    # pdb_id = "4GV1"  # 替换为实际的 PDB ID
    pdb_id  = sys.argv[1]
    name = sys.argv[2]
    save_dir = sys.argv[3]
    download_pdb_file(pdb_id, name, save_dir)

    # downloaded_file_path = download_pdb_file(pdb_id)
    # print(f"Downloaded file path: {downloaded_file_path}")
