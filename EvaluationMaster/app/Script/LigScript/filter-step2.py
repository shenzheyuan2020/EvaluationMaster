# filter_and_plot.py
# python filter.py /home/vesper/文档/lab/ligdown/p38a.csv p38a /home/vesper/文档/lab/ligdown/ 6.5

import sys
import pandas as pd
import matplotlib.pyplot as plt
import os

def filter_data(compounds_df, threshold):
    filtered_compounds = compounds_df[compounds_df['pIC50'] > threshold]
    return filtered_compounds

def plot_distribution(compounds_df, filter_name, save_dir, column='pIC50', title_suffix=''):
    plt.figure(figsize=(10, 8))
    plt.hist(compounds_df[column].dropna(), bins=50, color='skyblue', edgecolor='slategray', alpha=0.7)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(column)
    plt.ylabel('Frequency')
    title = f'{filter_name} Enhanced Histogram of {column} {title_suffix}'.strip()
    plt.title(title)
    
    pic50_to_ic50 = [
        (5, 10000),
        (6, 1000),
        (7, 100),
        (8, 10),
        (9, 1)  
    ]
    
    textstr = '\n'.join([f'pIC50 = {pic50}: IC50 = {ic50} nM' for pic50, ic50 in pic50_to_ic50])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right', bbox=props)
    
    image_file_path = os.path.join(save_dir, f"{filter_name}_filter.png")
    plt.savefig(image_file_path, bbox_inches='tight')
    plt.close()
    print(f"The enhanced histogram and comparison table have been saved to: {image_file_path}")


def main(csv_file_path, filter_name, save_dir, threshold):
    compounds_df = pd.read_csv(csv_file_path)
    filtered_compounds = filter_data(compounds_df, threshold)
    filtered_csv_file_path = os.path.join(save_dir, f"{filter_name}_filter.csv")  
    filtered_compounds.to_csv(filtered_csv_file_path, index=False)
    
    print(f"Filtered data saved to: {filtered_csv_file_path}")
    plot_distribution(filtered_compounds, filter_name, save_dir, column='pIC50', title_suffix='_filtered')


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python filter-step2.py <csv_file_path> <filter_name> <threshold> <save_dir>")
        sys.exit(1)

    csv_file_path = sys.argv[1]
    filter_name = sys.argv[2]
    threshold = float(sys.argv[3])  # Make sure this is a float since it's a threshold
    save_dir = sys.argv[4]          # This is the directory where the files will be saved
    main(csv_file_path, filter_name, save_dir, threshold)  # Arguments passed to main in correct order