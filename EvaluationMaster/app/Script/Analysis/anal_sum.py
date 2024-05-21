# import os
# import sys
# import pandas as pd

# def process_docking_results(base_path, file_type):
#     result_dfs = []

#     for root, _, files in os.walk(base_path):
#         for file in files:
#             if file.endswith('.xlsx') or file.endswith('.csv'):
#                 if 'docking_results.xlsx' in file:
#                     df = pd.read_excel(os.path.join(root, file))
#                 elif '_score.csv' in file:
#                     df = pd.read_csv(os.path.join(root, file))
#                 else:
#                     continue

#                 directory_name = os.path.basename(root)
#                 parts = directory_name.split('_')
#                 mol_name, software = parts[0], parts[1]

#                 if 'karma' in software.lower():
#                     score_columns = [col for col in df.columns if 'karma_score' in col and 'ff' not in col and 'aligned' not in col]
#                 else:
#                     score_columns = [col for col in df.columns if 'score' in col.lower()]

#                 rename_dict = {'Name': f'{mol_name}_{software}_Name'}
#                 for col in score_columns:
#                     rename_dict[col] = f'{mol_name}_{software}_Score'

#                 df = df.rename(columns=rename_dict)[[f'{mol_name}_{software}_Name'] + [f'{mol_name}_{software}_Score' for col in score_columns]]
#                 result_dfs.append(df)

#     if result_dfs:
#         combined_df = pd.concat(result_dfs, axis=1)
#     else:
#         combined_df = pd.DataFrame()

#     return combined_df

# def main():
#     if len(sys.argv) != 4:
#         print("Usage: script.py <decoy_dir> <inh_dir> <out_dir>")
#         sys.exit(1)

#     decoy_dir = sys.argv[1]
#     inh_dir = sys.argv[2]
#     out_dir = sys.argv[3]

#     decoy_results = process_docking_results(decoy_dir, 'decoy')
#     inh_results = process_docking_results(inh_dir, 'inh')

#     decoy_results.to_csv(os.path.join(out_dir, 'decoy_results.csv'), index=False)
#     inh_results.to_csv(os.path.join(out_dir, 'inh_results.csv'), index=False)

#     print("Decoy and inhibitor data have been processed and saved.")

# if __name__ == '__main__':
#     main()



import os
import sys
import pandas as pd

def process_docking_results(base_path, file_type):
    result_dfs = []

    for root, _, files in os.walk(base_path):
        for file in files:
            if file.endswith('.xlsx') or file.endswith('.csv'):
                if 'docking_results.xlsx' in file:
                    df = pd.read_excel(os.path.join(root, file))
                elif '_score.csv' in file:
                    df = pd.read_csv(os.path.join(root, file))
                else:
                    continue

                directory_name = os.path.basename(root)
                parts = directory_name.split('_')
                mol_name, software = parts[0], parts[1]

                if 'karma' in software.lower():
                    score_columns = [col for col in df.columns if 'karma_score' in col and 'ff' not in col and 'aligned' not in col]
                else:
                    score_columns = [col for col in df.columns if 'score' in col.lower()]

                rename_dict = {'Name': f'{mol_name}_{software}_Name'}
                for col in score_columns:
                    rename_dict[col] = f'{mol_name}_{software}_Score'

                df = df.rename(columns=rename_dict)
                
                # Ensure karma scores are negative
                if 'karma' in software.lower():
                    for col in score_columns:
                        if df[f'{mol_name}_{software}_Score'].gt(0).any():  # Check if there are any positive values
                            df[f'{mol_name}_{software}_Score'] = -df[f'{mol_name}_{software}_Score'].abs()

                # Select and order the columns to include in the result
                selected_columns = [f'{mol_name}_{software}_Name'] + [f'{mol_name}_{software}_Score' for col in score_columns]
                df = df[selected_columns]
                result_dfs.append(df)

    if result_dfs:
        combined_df = pd.concat(result_dfs, axis=1)
    else:
        combined_df = pd.DataFrame()

    return combined_df

def main():
    if len(sys.argv) != 4:
        print("Usage: script.py <decoy_dir> <inh_dir> <out_dir>")
        sys.exit(1)

    decoy_dir = sys.argv[1]
    inh_dir = sys.argv[2]
    out_dir = sys.argv[3]

    decoy_results = process_docking_results(decoy_dir, 'decoy')
    inh_results = process_docking_results(inh_dir, 'inh')

    decoy_results.to_csv(os.path.join(out_dir, 'decoy_results.csv'), index=False)
    inh_results.to_csv(os.path.join(out_dir, 'inh_results.csv'), index=False)

    print("Decoy and inhibitor data have been processed and saved.")

if __name__ == '__main__':
    main()
