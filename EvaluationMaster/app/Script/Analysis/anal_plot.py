import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from sklearn.metrics import roc_curve, auc
import seaborn as sns
import numpy as np
import os
import sys

def load_and_prepare_data(active_path, decoy_path):
    active_data = pd.read_csv(active_path)
    decoy_data = pd.read_csv(decoy_path)
    
    # Adjusting scores: Set all positive values to 0
    # Assuming the score columns end with '_Score'
    for column in active_data.columns:
        if '_Score' in column:
            active_data[column] = active_data[column].apply(lambda x: 0 if x > 0 else x)
    for column in decoy_data.columns:
        if '_Score' in column:
            decoy_data[column] = decoy_data[column].apply(lambda x: 0 if x > 0 else x)
    
    active_data['Label'] = 'Active'
    decoy_data['Label'] = 'Inactive'
    
    # Concatenate both datasets
    combined_data = pd.concat([active_data, decoy_data], ignore_index=True)
    return combined_data

def extract_names(merged_data):
    mol_names = [col.split('_')[0] for col in merged_data.columns if 'Score' in col]
    software_names = list(set(col.split('_')[1] for col in merged_data.columns if 'Score' in col))
    return mol_names, software_names

def perform_t_tests(merged_data, mol_names, software_names):
    ttest_results = {}
    for mol_name in set(mol_names):
        for software in software_names:
            score_col = f'{mol_name}_{software}_Score'
            name_col = f'{mol_name}_{software}_Name'
            # Select scores where 'Name' does not contain 'decoy'
            active_scores = merged_data[score_col][~merged_data[name_col].str.contains('decoy', na=False)]
            inactive_scores = merged_data[score_col][merged_data[name_col].str.contains('decoy', na=False)]
            if len(active_scores) > 0 and len(inactive_scores) > 0:
                t_stat, p_val = ttest_ind(active_scores.dropna(), inactive_scores.dropna())
                ttest_results[f'{mol_name}_{software}'] = (t_stat, p_val)
            else:
                ttest_results[f'{mol_name}_{software}'] = (np.nan, np.nan)
    return ttest_results


def save_ttest_results_to_csv(ttest_results, output_csv_path):
    # Convert the results to a DataFrame
    ttest_df = pd.DataFrame(ttest_results.items(), columns=['Combination', 'Statistics'])
    ttest_df[['T-Statistic', 'P-Value']] = pd.DataFrame(ttest_df['Statistics'].tolist(), index=ttest_df.index)
    ttest_df.drop('Statistics', axis=1, inplace=True)
    
    # Save the DataFrame to a CSV file
    ttest_df.to_csv(output_csv_path, index=False)
    print(f"Saved t-test results to {output_csv_path}")



# def load_and_prepare_data(active_path, decoy_path):
#     active_data = pd.read_csv(active_path)
#     decoy_data = pd.read_csv(decoy_path)
#     active_data['Label'] = 'Active'
#     decoy_data['Label'] = 'Inactive'
#     return pd.concat([active_data, decoy_data], ignore_index=True)

# def extract_names(merged_data):
#     mol_names = [col.split('_')[0] for col in merged_data.columns if 'Score' in col]
#     software_names = list(set(col.split('_')[1] for col in merged_data.columns if 'Score' in col))
#     return mol_names, software_names

# def perform_t_tests(merged_data, mol_names, software_names):
#     ttest_results = {}
#     for mol_name in set(mol_names):
#         for software in software_names:
#             score_col = f'{mol_name}_{software}_Score'
#             name_col = f'{mol_name}_{software}_Name'
#             try:
#                 active_scores = merged_data[score_col][merged_data[name_col].str.contains('CHEMBL', na=False)]
#                 inactive_scores = merged_data[score_col][merged_data[name_col].str.contains('decoy', na=False)]
#                 if len(active_scores) > 0 and len(inactive_scores) > 0:
#                     t_stat, p_val = ttest_ind(active_scores.dropna(), inactive_scores.dropna())
#                     ttest_results[f'{mol_name}_{software}'] = (t_stat, p_val)
#                 else:
#                     ttest_results[f'{mol_name}_{software}'] = (np.nan, np.nan)
#                     print(f"No data for t-test between active and inactive molecules for {mol_name} with {software}.")
#             except Exception as e:
#                 ttest_results[f'{mol_name}_{software}'] = (np.nan, np.nan)
#                 print(f"Error performing t-test for {mol_name} with {software}: {e}")
#     return ttest_results

# def save_ttest_results_to_csv(ttest_results, output_csv_path):
#     ttest_df = pd.DataFrame([(k, *v) for k, v in ttest_results.items()], columns=['Combination', 'T-Statistic', 'P-Value'])
#     ttest_df.to_csv(output_csv_path, index=False)
#     print(f"Saved t-test results to {output_csv_path}")









def plot_roc_curves(merged_data, mol_names, software_names, output_path):
    plt.figure(figsize=(10, 8))
    for software in software_names:
        for mol_name in set(mol_names):
            score_col = f'{mol_name}_{software}_Score'
            if score_col in merged_data.columns:
                scores = merged_data[score_col].replace([np.inf, -np.inf], np.nan).dropna()
                labels = merged_data['Label'].loc[scores.index].apply(lambda x: 0  if x == 'Active' else 1)
                fpr, tpr, _ = roc_curve(labels, scores)
                roc_auc = auc(fpr, tpr)
                plt.plot(fpr, tpr, label=f'{mol_name}_{software} (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], 'k--', label='Random Classifier (AUC = 0.50)')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for Docking Software')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'roc_plot.png'),dpi=300)
    plt.close()


# def plot_boxplots(merged_data, mol_names, software_names, output_path):
#     num_softwares = len(software_names)
#     fig, axes = plt.subplots(1, num_softwares, figsize=(15, 10), squeeze=False)
#     axes = axes.flatten()
#     for idx, software in enumerate(sorted(software_names)):
#         ax = axes[idx]
#         for mol_name in sorted(mol_names):
#             score_col = f'{mol_name}_{software}_Score'
#             if score_col in merged_data.columns:
#                 sns.boxplot(x='Label', y=score_col, hue='Label', data=merged_data, ax=ax, palette='pastel')
#                 ax.set_title(f'Boxplot for {software}')
#     plt.tight_layout()
#     plt.savefig(os.path.join(output_path, 'boxplot.png'), dpi=300)
#     plt.close()


def plot_boxplots(merged_data, mol_names, software_names, output_path):
    total_plots = len(mol_names) * len(software_names)
    # 计算接近平方根的行列数
    num_rows = int(np.ceil(np.sqrt(total_plots)))
    num_cols = int(np.ceil(total_plots / num_rows))
    # num_rows = 7
    # num_cols = 3
    # 创建子图
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows), squeeze=False)
    plot_idx = 0
    for software in sorted(software_names):
        for mol_name in sorted(mol_names):
            row, col = divmod(plot_idx, num_cols)
            ax = axes[row, col]
            score_col = f'{mol_name}_{software}_Score'
            if score_col in merged_data.columns:
                sns.boxplot(x='Label', y=score_col, hue='Label', data=merged_data, ax=ax, palette='pastel')
                ax.set_title(f'Boxplot for {mol_name} - {software}', fontsize=16, fontweight='bold')  # 调整标题字体大小和加粗
                ax.set_xlabel('Label', fontsize=16, fontweight='bold')  # 调整X轴标签大小和加粗
                ax.set_ylabel(score_col, fontsize=16, fontweight='bold')  # 调整Y轴标签大小和加粗
                ax.tick_params(axis='both', labelsize=9)  # 调整轴标记的字体大小
            else:
                ax.set_title(f'No data for {mol_name} - {software}', fontsize=16, fontweight='bold')  # 同样调整无数据时的标题
                ax.set_visible(False)  # 如果没有数据，隐藏这个子图
            plot_idx += 1
    
    # 用于隐藏多余的子图
    for i in range(plot_idx, num_rows * num_cols):
        row, col = divmod(i, num_cols)
        axes[row, col].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'boxplot.png'), dpi=300)
    plt.close()


import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

def plot_violins(merged_data, mol_names, software_names, output_path):
    total_plots = len(mol_names) * len(software_names)
    # 计算接近平方根的行列数
    num_rows = int(np.ceil(np.sqrt(total_plots)))
    num_cols = int(np.ceil(total_plots / num_rows))
    # num_rows = 7
    # num_cols = 3
    # 创建子图
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows), squeeze=False)
    plot_idx = 0
    for software in sorted(software_names):
        for mol_name in sorted(mol_names):
            row, col = divmod(plot_idx, num_cols)
            ax = axes[row, col]
            score_col = f'{mol_name}_{software}_Score'
            if score_col in merged_data.columns:
                sns.violinplot(x='Label', y=score_col, hue='Label', data=merged_data, ax=ax, palette='pastel')
                ax.set_title(f'Violin plot for {mol_name} - {software}', fontsize=12, fontweight='bold')
                ax.set_xlabel('Label', fontsize=10, fontweight='bold')
                ax.set_ylabel(score_col, fontsize=10, fontweight='bold')
                ax.tick_params(axis='both', labelsize=9)
                ax.legend(title='Label', title_fontsize='10', fontsize='9', frameon=False)  # Customize legend
            else:
                ax.set_visible(False)  # 如果没有数据，隐藏这个子图
            plot_idx += 1
    
    # 用于隐藏多余的子图
    for i in range(plot_idx, num_rows * num_cols):
        row, col = divmod(i, num_cols)
        axes[row, col].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'violinplot.png'), dpi=300)
    plt.close()

# Example of function call
# plot_violins(merged_data, ['Molecule1', 'Molecule2'], ['Software1', 'Software2'], 'path_to_output')


def plot_box_individual(merged_data, mol_names, software_names, output_path):
    for software in sorted(software_names):
        plt.figure(figsize=(6, 4))  
        for mol_name in mol_names:
            score_col = f'{mol_name}_{software}_Score'
            if score_col in merged_data.columns:
                sns.boxplot(x='Label', y=score_col, data=merged_data, palette='pastel')
                plt.title(f'box plot for {software}')
                plt.ylabel('Score')
                plt.xlabel('Label')
                plt.tight_layout()
                plt.savefig(os.path.join(output_path, f'boxplot_{software}_{mol_name}.png'), dpi=300)
                plt.close()

def plot_violins_individual(merged_data, mol_names, software_names, output_path):
    for software in sorted(software_names):
        plt.figure(figsize=(6, 4))  
        for mol_name in mol_names:
            score_col = f'{mol_name}_{software}_Score'
            if score_col in merged_data.columns:
                sns.violinplot(x='Label', y=score_col, data=merged_data, palette='pastel')
                plt.title(f'Violin Plot for {software}')
                plt.ylabel('Score')
                plt.xlabel('Label')
                plt.tight_layout()
                plt.savefig(os.path.join(output_path, f'violinplot_{software}_{mol_name}.png'), dpi=300)
                plt.close()
 
# def plot_binned_bar_charts(merged_data, mol_names, software_names, output_path):
#     for software in software_names:
#         for mol_name in set(mol_names):
#             score_col = f'{mol_name}_{software}_Score'
#             if score_col in merged_data.columns:
#                 # Binning the score data
#                 merged_data['Binned'] = (merged_data[score_col] // 0.5) * 0.5
#                 grouped = merged_data.groupby(['Binned', 'Label']).size().unstack().fillna(0)
                
#                 # Normalizing the data
#                 total_label_active = grouped['Active'].sum()
#                 total_label_inactive = grouped['Inactive'].sum()
#                 grouped['Active'] = grouped['Active'] / total_label_active
#                 grouped['Inactive'] = grouped['Inactive'] / total_label_inactive
                
#                 # Plotting the histograms
#                 fig, ax = plt.subplots(figsize=(15,7))
#                 x = np.arange(len(grouped.index))
#                 width = 0.35  # the width of the bars
                
#                 ax.bar(x - width/2, grouped['Active'], width, label='Active', color='blue')
#                 ax.bar(x + width/2, grouped['Inactive'], width, label='Inactive', color='red')
                
#                 # Add some text for labels, title, and custom x-axis tick labels
#                 ax.set_ylabel('Fraction of Molecules')
#                 ax.set_xlabel('Binned Score Value')
#                 ax.set_title(f'Fraction of Molecules by Binned Score Value for {software}')
#                 ax.set_xticks(x)
#                 ax.set_xticklabels(grouped.index, rotation=45)
#                 ax.legend()
                
#                 # Save the plot
#                 plt.tight_layout()
#                 plt.savefig(os.path.join(output_path, f'frequency_histogram_{software}.png'), dpi=300)
#                 plt.close()


def plot_binned_bar_charts(merged_data, mol_names, software_names, output_path):
    total_plots = len(mol_names) * len(software_names)
    num_rows = int(np.ceil(np.sqrt(total_plots)))
    num_cols = int(np.ceil(total_plots / num_rows))
    # num_rows = 7
    # num_cols = 3
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows), squeeze=False)
    plot_idx = 0
    
    colors = {'Active': '#add8e6', 'Inactive': '#ffcc99'}  # Define colors for the bars
    
    for software in sorted(software_names):
        for mol_name in sorted(set(mol_names)):
            row, col = divmod(plot_idx, num_cols)
            ax = axes[row, col]
            score_col = f'{mol_name}_{software}_Score'
            
            if score_col in merged_data.columns:
                # Process scores: set all positive scores to zero
                modified_scores = merged_data[score_col].apply(lambda x: 0 if x > 0 else x)
                min_score = int(np.floor(modified_scores.min()))
                max_score = 0  # Max score set to zero since we're changing positive values to zero

                # Define bins, ensuring we include zero if min_score is negative
                bins = np.linspace(min_score, max_score, 21) if min_score < 0 else np.array([min_score])

                merged_data['Binned'] = pd.cut(modified_scores, bins=bins, labels=bins[:-1], right=False)
                grouped = merged_data.groupby(['Binned', 'Label']).size().unstack().fillna(0)
                
                # Normalize data
                total_label_active = grouped['Active'].sum()
                total_label_inactive = grouped['Inactive'].sum()
                grouped['Active'] = grouped['Active'] / total_label_active
                grouped['Inactive'] = grouped['Inactive'] / total_label_inactive
                
                x = np.arange(len(grouped.index))  # Index for bars
                width = 0.35  # Bar width
                ax.bar(x - width/2, grouped['Active'], width, label='Active', color=colors['Active'])
                ax.bar(x + width/2, grouped['Inactive'], width, label='Inactive', color=colors['Inactive'])
                
                ax.set_ylabel('Fraction of Molecules')
                ax.set_xlabel('Binned Score Value')
                ax.set_title(f'{mol_name} - {software}')
                ax.set_xticks(x)
                ax.set_xticklabels([f'{b:.1f}' for b in bins[:-1]], rotation=45)
                ax.legend()
            else:
                ax.set_visible(False)  # Hide axis if no data
            plot_idx += 1
    
    # Hide unused subplots
    for i in range(plot_idx, num_rows * num_cols):
        row, col = divmod(i, num_cols)
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'binned_bar_charts.png'), dpi=300)
    plt.close()


def main():
    if len(sys.argv) < 4:
        print("Usage: python script.py <active_csv_path> <decoy_csv_path> <output_directory>")
        return
    decoy_path = sys.argv[1]
    active_path = sys.argv[2]
    output_path = sys.argv[3]

    merged_data = load_and_prepare_data(active_path, decoy_path)
    mol_names, software_names = extract_names(merged_data)
    ttest_results = perform_t_tests(merged_data, mol_names, software_names)
    for key, value in ttest_results.items():
        print(f'{key}: t-statistic = {value[0]}, p-value = {value[1]}')
    output_csv_path = os.path.join(output_path, 'ttest_results.csv')
    save_ttest_results_to_csv(ttest_results, output_csv_path)
    
    plot_roc_curves(merged_data, mol_names, software_names, output_path)
    plot_boxplots(merged_data, mol_names, software_names, output_path)
    plot_violins(merged_data, mol_names, software_names, output_path)
    plot_binned_bar_charts(merged_data, mol_names, software_names, output_path)
    plot_violins_individual(merged_data, mol_names, software_names, output_path)
    plot_box_individual(merged_data, mol_names, software_names, output_path)
if __name__ == '__main__':
    main()

