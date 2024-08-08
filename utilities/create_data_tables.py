import os
import pandas as pd
import numpy as np
from openpyxl import load_workbook
from sklearn.preprocessing import StandardScaler

def process_reg_mat(reg_str, phospho, targets, cluster_mat, TF, candidates):
    pass

def main():
    current_directory = os.getcwd()

    target_str = "JA_FPKM.csv"
    cluster_str = "JA_FPKM_means.csv"
    TF_str = "Arabidopsis_TFs_AGRIS.xlsx"
    candidates_str = "SC-ION_candidates.txt"
    
    # Read in files
    target_mat = pd.read_csv(f"{current_directory}\\test_files\\Raw_Files\{target_str}", index_col=0)
    cluster_mat = pd.read_csv(f"{current_directory}\\test_files\\Raw_Files\{cluster_str}", index_col=0)
    cluster_mat = cluster_mat[sorted(cluster_mat.columns)]
    TF = pd.read_excel(f"{current_directory}\\test_files\\Raw_Files\{TF_str}")
    candidates = pd.read_csv(f"{current_directory}\\test_files\\Raw_Files\{candidates_str}", header=None)

    # Process target_mat first
    # Here, we only take the 0 and 2hr timepoints
    my_target_mat = target_mat.filter(regex='hr0_|hr2_')
    scaler = StandardScaler()
    my_target_mat = pd.DataFrame(scaler.fit_transform(my_target_mat.T).T, index=my_target_mat.index, columns=my_target_mat.columns)
    # Remove rows where all values are 0
    my_target_mat = my_target_mat.loc[~(my_target_mat == 0.0).all(axis=1)]
    my_target_mat.dropna(inplace=True)
    my_target_mat.to_csv('target_mat_RNA.csv', na_rep='0')

    # Save target list
    my_targets = my_target_mat.index[my_target_mat.index.isin(candidates[0])]
    my_targets_df = pd.DataFrame(my_targets, columns=['Targets'])
    my_targets_df.to_csv('target_list_RNA.csv', index=False)

if __name__ == "__main__":
    main()