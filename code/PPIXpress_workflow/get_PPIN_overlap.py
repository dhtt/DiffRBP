#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
from functools import reduce
from matplotlib import pyplot as plt

def jaccard_similarity(list1, list2):
    list1 = set(list1) - set(['', 'nan'])
    list2 = set(list2) - set(['', 'nan'])
    if len(list1.union(list2)) == 0:
        return 0
    else :
        return(len(list1.intersection(list2)) / len(list1.union(list2)))


def get_PPIN_overlap(PPIXpress_path, result_path):
    sim_plot_path = os.path.join(result_path, 'PPIN_sim_cluster.png')

    PPIN_files = []
    for root, dirs, files in os.walk(PPIXpress_path):
        for file in files:
            if file.endswith('_ppin.txt') & (not file.startswith('.')) & (not file.startswith('reference')):
                PPIN_files.append(os.path.join(root, file))
    PPIN_files = sorted(PPIN_files, key=str.lower)
    
    PPIN = dict()
    prev_tissue = ''
    idx = 1
    for file in PPIN_files:
        current_tissue = file.split('/')[-2]
        idx += 1 if current_tissue != prev_tissue else 1 
        
        df = pd.read_csv(file, sep=' ', header=0)
        PPIN[f"{current_tissue}_{idx}"] = df['Protein1'].str.cat(df['Protein2'], sep='_').to_list()
        prev_tissue = current_tissue
    
    PPIN_similarity = [[x, y, jaccard_similarity(PPIN[x], PPIN[y])] for x in PPIN.keys() for y in PPIN.keys()]
    PPIN_similarity = pd.DataFrame(PPIN_similarity, columns=['PPIN1', 'PPIN2', 'Jaccard_similarity']).pivot(index='PPIN1', columns='PPIN2', values='Jaccard_similarity')
    
    sim_plot = sns.clustermap(PPIN_similarity, linewidths=0.5, figsize=(18, 16), cmap='viridis')
    sim_plot.savefig(sim_plot_path, dpi=300)

                     
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--PPIXpress_path', dest='PPIXpress_path', type=str, help='Path to all results PPIXpress analysis.')
    parser.add_argument('--result_path', dest='result_path', type=str, help='Path to all results of PPIN_overlap analysis.')
    args = parser.parse_args()

    PPIXpress_path = args.PPIXpress_path
    result_path = args.result_path
    
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    get_PPIN_overlap(PPIXpress_path=PPIXpress_path, result_path=result_path)