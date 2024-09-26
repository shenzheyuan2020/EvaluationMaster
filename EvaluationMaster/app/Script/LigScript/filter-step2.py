# filter_and_plot.py
# python filter.py /home/vesper/文档/lab/ligdown/p38a.csv p38a /home/vesper/文档/lab/ligdown/ 6.5

# import sys
# import pandas as pd
# import matplotlib.pyplot as plt
# import os

# def filter_data(compounds_df, threshold):
#     filtered_compounds = compounds_df[compounds_df['pIC50'] > threshold]
#     return filtered_compounds

# def plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix=''):
#     plt.figure(figsize=(10, 8))
#     plt.hist(compounds_df[column].dropna(), bins=50, color='skyblue', edgecolor='slategray', alpha=0.7)
#     plt.grid(axis='y', alpha=0.75)
#     plt.xlabel(column)
#     plt.ylabel('Frequency')
#     title = f'{filter_name} Enhanced Histogram of {column} {title_suffix}'.strip()
#     plt.title(title)
    
#     pic50_to_ic50 = [
#         (5, 10000),
#         (6, 1000),
#         (7, 100),
#         (8, 10),
#         (9, 1)  
#     ]
    
#     textstr = '\n'.join([f'pIC50 = {pic50}: IC50 = {ic50} nM' for pic50, ic50 in pic50_to_ic50])
#     props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#     plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
#              verticalalignment='top', horizontalalignment='right', bbox=props)
    
#     image_file_path = os.path.join(save_dir, f"{filter_name}_filter.png")
#     plt.savefig(image_file_path, bbox_inches='tight')
#     plt.close()
#     print(f"The enhanced histogram and comparison table have been saved to: {image_file_path}")


# def main(csv_file_path, filter_name, save_dir, threshold):
#     compounds_df = pd.read_csv(csv_file_path)
#     filtered_compounds = filter_data(compounds_df, threshold)
#     filtered_csv_file_path = os.path.join(save_dir, f"{filter_name}_filter.csv")  
#     filtered_compounds.to_csv(filtered_csv_file_path, index=False)
    
#     print(f"Filtered data saved to: {filtered_csv_file_path}")
#     plot_distribution(filtered_compounds, filter_name, save_dir, column='pIC50', title_suffix='_filtered')


# if __name__ == "__main__":
#     if len(sys.argv) != 5:
#         print("Usage: python filter-step2.py <csv_file_path> <filter_name> <threshold> <save_dir>")
#         sys.exit(1)

#     csv_file_path = sys.argv[1]
#     filter_name = sys.argv[2]
#     threshold = float(sys.argv[3])  # Make sure this is a float since it's a threshold
#     save_dir = sys.argv[4]          # This is the directory where the files will be saved
#     main(csv_file_path, filter_name, save_dir, threshold)  # Arguments passed to main in correct order



import sys
import pandas as pd
import matplotlib.pyplot as plt
import os

def filter_data(compounds_df, selected_type, threshold):
    """根据给定阈值过滤化合物数据。"""
    filtered_compounds = compounds_df.copy()
    if selected_type in filtered_compounds.columns:
        filtered_compounds = filtered_compounds[filtered_compounds[selected_type] > threshold]
    return filtered_compounds

def plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix=''):
    """绘制指定列的分布直方图。"""
    plt.figure(figsize=(10, 8))
    plt.hist(compounds_df[column].dropna(), bins=50, color='skyblue', edgecolor='slategray', alpha=0.7)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(column)
    plt.ylabel('Frequency')
    title = f'{filter_name} Enhanced Histogram of {column} {title_suffix}'.strip()
    plt.title(title)

    # 根据不同的列动态生成对照表
    if column == 'pIC50':
        pic50_to_ic50 = [(5, 10000), (6, 1000), (7, 100), (8, 10), (9, 1)]
        textstr = '\n'.join([f'pIC50 = {pic50}: IC50 = {ic50} nM' for pic50, ic50 in pic50_to_ic50])
    elif column == 'pKd':
        kd_to_pkd = [(5, 10000), (6, 1000), (7, 100), (8, 10), (9, 1)]
        textstr = '\n'.join([f'pKd = {pkd}: Kd = {kd} nM' for pkd, kd in kd_to_pkd])
    elif column == 'pKi':
        ki_to_pki = [(5, 10000), (6, 1000), (7, 100), (8, 10), (9, 1)]
        textstr = '\n'.join([f'pKi = {pki}: Ki = {ki} nM' for pki, ki in ki_to_pki])

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right', bbox=props)

    image_file_path = os.path.join(save_dir, f"{filter_name}_{column}_filter.png")
    plt.savefig(image_file_path, bbox_inches='tight')
    plt.close()
    print(f"The enhanced histogram and comparison table have been saved to: {image_file_path}")

def main(csv_file_path, filter_name, save_dir, selected_type, threshold):
    """主函数：读取数据，过滤并绘图。"""
    compounds_df = pd.read_csv(csv_file_path)
    filtered_compounds = filter_data(compounds_df, selected_type, threshold)
    

    filtered_csv_file_path = os.path.join(save_dir, f"{filter_name}_{selected_type}_filter.csv")  
    filtered_compounds.to_csv(filtered_csv_file_path, index=False)
    
    print(f"Filtered data saved to: {filtered_csv_file_path}")

    # 绘制过滤后的列的直方图
    if selected_type in filtered_compounds.columns:
        plot_distribution(filtered_compounds, filter_name, save_dir, column=selected_type, title_suffix='_filtered')

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python filter.py <csv_file_path> <filter_name> <save_dir> <column> <threshold>")
        print("Example: python filter.py /path/to/p38a.csv p38a /path/to/save/ pIC50 6.0")
        sys.exit(1)

    csv_file_path = sys.argv[1]
    filter_name = sys.argv[2]
    save_dir = sys.argv[3]
    selected_type = sys.argv[4]  # 需要过滤的列名
    threshold = float(sys.argv[5])  # 过滤阈值

    main(csv_file_path, filter_name, save_dir, selected_type, threshold)
