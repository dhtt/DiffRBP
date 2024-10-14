#!/usr/bin/env python3

import pandas as pd 
from os import path, makedirs, listdir
from get_nonPPIN import *
from annotate_by_mRNPchrono import *
from analyze_distribution import *
import argparse

def main(result_path: str, PPIN_path: str, path_to_mRNPchrono_file: str):
    makedirs(result_path, exist_ok=True)
    tissue = PPIN_path.split('/')[-1]
    sample_ids = [x.split('/')[-1].split('_')[0] for x in listdir(PPIN_path) if x.endswith('_ppin.txt') & (not x.startswith('.'))]

    for sample_no in sample_ids:    
        path_to_sample_PPIN = f"{PPIN_path}/{sample_no}_ppin.txt"
            
        # Step 1: Annotate PPIs with mRNPchrono
        PPIN = pd.read_csv(path_to_sample_PPIN, sep=' ', header=0)
        mRNPchrono_table = pd.read_excel(path_to_mRNPchrono_file, header=0)
        mRNPchrono_table = get_compartment_list(mRNPchrono_table)

        PPIN_df = annotate_by_mRNPchrono(PPIN=PPIN, mRNPchrono_table=mRNPchrono_table, cleanup=True)
        PPIN_df['is_PPI'] = 1
        # print(PPIN_df.head())
        
        # Step 2: Get non-PPI set
        nonPPIN_df = get_nonPPIN(PPIN=PPIN_df, mRNPchrono_table=mRNPchrono_table, cleanup=True)
        nonPPIN_df['is_PPI'] = 0
        # print(nonPPIN_df.head())

        annotated_PPIN = pd.concat([PPIN_df, nonPPIN_df])
        annotated_PPIN['tissue'] = tissue
        annotated_PPIN['sample'] = sample_no
        
        # Step 3: Prepare data for analysis of binding time distribution
        annotated_PPIN = get_cluster_distance(annotated_PPIN)
        annotated_PPIN = get_cluster_overlap(annotated_PPIN)
        annotated_PPIN = get_binding_time_difference(annotated_PPIN)
        annotated_PPIN = get_compartment_overlap(annotated_PPIN)

        annotated_PPIN.to_csv(path.join(result_path, f"annotated_PPIN_{sample_no}.csv"), index=False, sep='\t')
        # print(annotated_PPIN.head())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--PPIN_path', dest='PPIN_path', type=str, help='Path to results of PPIXpress analysis.')
    parser.add_argument('--result_path', dest='result_path', type=str, help='Path to all results of PPIN_vs_binding_time analysis.')
    parser.add_argument('--mRNP_chronology_table_path', dest='mRNP_chronology_table_path', type=str, help='Path to mRNP_chronology_hsa_HeLa-v1.xlsx.')
    args = parser.parse_args()

    # tissue = args.tissue
    # sample_no =  args.sample_no
    # PPIN_path = f"../../data/RBP2GO/PPIXpress/RBP2GO_all_options/ResultFiles_renamed/"
    # result_path = f"../../data/RBP2GO/PPIN_vs_binding_time/{tissue}/"
    # path_to_sample_PPIN = f"{args.PPIN_path}/{tissue}/{sample_no}_ppin.txt"
    PPIN_path = args.PPIN_path
    result_path = args.result_path
    mRNP_chronology_table_path = args.mRNP_chronology_table_path

    main(PPIN_path=PPIN_path, path_to_mRNPchrono_file=mRNP_chronology_table_path, result_path=result_path)