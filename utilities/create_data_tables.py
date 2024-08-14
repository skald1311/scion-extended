import os
import pandas as pd
import numpy as np
from openpyxl import load_workbook
from sklearn.preprocessing import StandardScaler

# Function to process regulator matrix
def process_reg_mat(reg_str, phospho, targets, cluster_mat, TF, candidates):
    reg_mat = pd.read_csv(reg_str, index_col=0)
    scaler = StandardScaler()
    if phospho:
        # Merge multiplicities
        reg_mat_grouped = reg_mat.groupby('original.id').mean(numeric_only=True)
        my_data = reg_mat_grouped.filter(regex='JA|original.id')
        
        # Get site information for each ID
        reg_mat['original.id'] = reg_mat['original.id'].astype(float)
        reg_mat.sort_values(by='original.id', inplace=True)
        
        my_ids = reg_mat['original.id'].drop_duplicates().values
        print(reg_mat['Protein'].values)
        print(len(set(reg_mat['Protein'].values)))
        my_proteins = reg_mat['Protein'].drop_duplicates().values
        my_sites = reg_mat['Positions.within.proteins'].drop_duplicates().values
        #print(len(my_ids))
        #print(len(my_proteins))
        #print(len(my_sites))
        my_names = []
        for i in range(len(my_proteins)):
            my_pro = my_proteins[i].split('.')[0]
            my_site = my_sites[i].split(';')[0]
            my_id = f"{my_pro}.p{my_site}"
            my_names.append([my_id, my_ids[i]])
        
        my_names_df = pd.DataFrame(my_names, columns=['Protein.Site', 'original.id'])
        
        final_data = pd.merge(my_names_df, my_data, on='original.id', sort=False)
        
        # # Collapse duplicate sites
        # if final_data.duplicated('Protein.Site').any():
        #     final_data_grouped = final_data.groupby('Protein.Site').mean()
        #     my_pros = final_data_grouped.index
        #     final_data = final_data_grouped.drop(columns=['original.id'])
        #     final_data.index = my_pros
        # else:
        #     final_data.index = final_data['Majority.protein.IDs']
        #     final_data.drop(columns=['original.id', 'Majority.protein.IDs'], inplace=True)
        
        # final_data = pd.DataFrame(scaler.fit_transform(final_data), index=final_data.index, columns=final_data.columns)
        # final_data.to_csv('reg_mat_phospho.csv', na_rep='0')
        
        # TF_list = TF[TF.iloc[:, 0].isin(candidates.iloc[:, 0])]
        # my_genes = [gene for genes in final_data.index.str.split('.') for gene in genes if 'AT' in gene]
        # TF_final = final_data.index[final_data.index.str.split('.').map(set(my_genes).intersection)]
        
        # pd.DataFrame(TF_final).to_csv('reg_list_phospho.csv', index=False, header=False)
    
        

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


    # # Process proteome data
    # proteome_str = "Normalized_values_protein.csv"
    # process_reg_mat(f"{current_directory}\\test_files\\Raw_Files\{proteome_str}", False, my_targets, cluster_mat, TF, candidates)

    # Process phospho data
    phospho_str = "Normalized_values_phospho.csv"
    process_reg_mat(f"{current_directory}\\test_files\\Raw_Files\{phospho_str}", True, my_targets, cluster_mat, TF, candidates)

if __name__ == "__main__":
    main()