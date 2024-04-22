import pandas as pd
from rdkit import Chem
import os

import pandas as pd
from rdkit import Chem
import os

def is_smiles(smiles_string):
    """
    Check if a string is a valid smiles using RDKit.
    """
    try:
        mol = Chem.MolFromsmiles(smiles_string)
        return mol is not None
    except:
        return False

def detect_delimiter(file_path):
    """
    Detects the delimiter in a given file (tab, comma, or space).
    """
    with open(file_path, 'r') as file:
        line = file.readline()
        if '\t' in line:
            return '\t'
        elif ',' in line:
            return ','
        else:
            return ' '

def process_smi_file(smi_file):
    """
    Processes a .smi file and converts it to a standard CSV format.
    """
    delimiter = detect_delimiter(smi_file)
    df = pd.read_csv(smi_file, delimiter=delimiter, header=None, nrows=2)
    num_columns = len(df.columns)

    if num_columns == 1:
        # Single column .smi file
        df = pd.read_csv(smi_file, delimiter=delimiter, header=None)
        df.columns = ['smiles']
    elif num_columns >= 2:
        # File with two or more columns; using the second column as smiles
        df = pd.read_csv(smi_file, delimiter=delimiter, header=None)
        df.columns = ['Extra', 'smiles'] + [f'Extra_{i}' for i in range(2, num_columns)]
        df = df[['smiles']]
    else:
        raise ValueError("Unsupported number of columns in .smi file.")

    df['Name'] = ['lig_' + str(i) for i in range(1, len(df) + 1)]

    standard_csv_file = os.path.splitext(smi_file)[0] + '_standard.csv'
    df.to_csv(standard_csv_file, index=False)
    return standard_csv_file



def process_csv_file(csv_file):
    """
    Processes a CSV file and converts it to a standard CSV format with smiles and Name columns.
    """
    df = pd.read_csv(csv_file)

    if 'smiles' not in df.columns:
        raise ValueError("CSV file must contain a 'smiles' column.")

    new_df = pd.DataFrame()
    new_df['smiles'] = df['smiles']
    new_df['Name'] = ['lig_' + str(i) for i in range(1, len(df) + 1)]

    standard_csv_file = os.path.splitext(csv_file)[0] + '_standard.csv'
    new_df.to_csv(standard_csv_file, index=False)
    return standard_csv_file

def convert_file_to_csv(input_file):
    """
    Converts various file formats to a standard CSV format.
    """
    file_ext = os.path.splitext(input_file)[1].lower()

    if file_ext in ['.xlsx', '.xls']:
        df = pd.read_excel(input_file)
        standard_csv_file = os.path.splitext(input_file)[0] + '_standard.csv'
        df['Name'] = ['lig_' + str(i) for i in range(1, len(df) + 1)]
        df[['smiles', 'Name']].to_csv(standard_csv_file, index=False)
    elif file_ext == '.txt':
        df = pd.read_csv(input_file)
        standard_csv_file = os.path.splitext(input_file)[0] + '_standard.csv'
        df['Name'] = ['lig_' + str(i) for i in range(1, len(df) + 1)]
        df[['smiles', 'Name']].to_csv(standard_csv_file, index=False)
    elif file_ext == '.csv':
        standard_csv_file = process_csv_file(input_file)
    elif file_ext == '.smi':
        standard_csv_file = process_smi_file(input_file)
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    return standard_csv_file

# Example usage:
# input_file = 'path/to/your/input/file.smi'  # Replace with actual file path
# csv_file = convert_file_to_csv(input_file)
# print(f"Converted to standard CSV: {csv_file}")
