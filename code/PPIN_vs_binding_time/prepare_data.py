import pandas as pd 
from os import path, makedirs
from get_nonPPIN import *
from annotate_by_mRNPchrono import *
from analyze_distribution import *
import seaborn as sns
import matplotlib.pyplot as plt


def main(path_to_PPI_file: str, path_to_mRNPchrono_file: str, result_path: str):
    # Step 1: Annotate PPIs with mRNPchrono
    PPIN = pd.read_csv(path_to_PPI_file, sep=' ', header=0)
    mRNPchrono_table = pd.read_excel(path_to_mRNPchrono_file, header=0)
    mRNPchrono_table = get_compartment_list(mRNPchrono_table)

    PPIN_df = annotate_by_mRNPchrono(PPIN=PPIN, mRNPchrono_table=mRNPchrono_table, cleanup=True)
    PPIN_df['is_PPI'] = 1
    print(PPIN_df.head())
    
    # Step 2: Get non-PPI set
    nonPPIN_df = get_nonPPIN(PPIN=PPIN_df, mRNPchrono_table=mRNPchrono_table, cleanup=True)
    nonPPIN_df['is_PPI'] = 0
    print(nonPPIN_df.head())

    annotated_PPIN = pd.concat([PPIN_df, nonPPIN_df])
    annotated_PPIN['tissue'] = f'{tissue}'
    annotated_PPIN['sample'] = f'{sample_no}'
    
    # Step 3: Prepare data for analysis of binding time distribution
    annotated_PPIN = get_cluster_distance(annotated_PPIN)
    annotated_PPIN = get_cluster_overlap(annotated_PPIN)
    annotated_PPIN = get_binding_time_difference(annotated_PPIN)
    annotated_PPIN = get_compartment_overlap(annotated_PPIN)

    annotated_PPIN.to_csv(path.join(result_path, f"annotated_PPIN_{sample_no}.csv"), index=False, sep='\t')
    print(annotated_PPIN.head())


if __name__ == '__main__':
    tissue = 'adipose'
    sample_no = 2
    path_to_RBP2GO_PPINs = f"../../data/RBP2GO/PPIXpress/RBP2GO_all_options/ResultFiles_renamed/{tissue}/"
    path_to_sample_PPIN = path_to_RBP2GO_PPINs + f"{sample_no}_ppin.txt"
    result_path = f"../../data/RBP2GO/PPIN_vs_binding_time/{tissue}/"
    makedirs(result_path, exist_ok=True)
    
    main(path_to_PPI_file=path_to_sample_PPIN, 
         path_to_mRNPchrono_file="../../data/mRNP_chrono/mRNP_chronology_hsa_HeLa-v1.xlsx",
         result_path=result_path)