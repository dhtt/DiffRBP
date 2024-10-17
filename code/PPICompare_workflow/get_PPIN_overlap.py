#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
from functools import reduce
from matplotlib import pyplot as plt
import numpy as np

def jaccard_similarity(list1, list2):
    list1 = set(list1) - set(['', 'nan'])
    list2 = set(list2) - set(['', 'nan'])
    if len(list1.union(list2)) == 0:
        return 0
    else :
        return(len(list1.intersection(list2)) / len(list1.union(list2)))

def get_PPIN_overlap(PPICompare_path, result_path):
    sim_plot_path = os.path.join(result_path, 'PPIN_sim_cluster.png')
    diff_sim_plot_path = os.path.join(result_path, 'diff_PPIN_sim_cluster.png')

    PPIN_files = []
    for root, dirs, files in os.walk(PPICompare_path):
        for file in files:
            if file == 'differential_network.txt':
                PPIN_files.append(os.path.join(root, file))
    PPIN_files = sorted(PPIN_files, key=str.lower)
    
    PPINs, diff_PPINs = [], dict()
    for file in PPIN_files:
        tissue_1 = file.split('/')[-2].split('_')[0]
        tissue_2 = file.split('/')[-2].split('_')[1]
        df = pd.read_csv(file, sep=' ', header=0)
        df['p-val_adj'] = df['p-val_adj'].astype(float)
        df = df[df['p-val_adj'] < 0.05]
        
        if df.shape[0] > 0:
            # Pairwise tissue PPIN
            PPINs.append([tissue_1, tissue_2, df.shape[0]])
            PPINs.append([tissue_2, tissue_1, df.shape[0]])
                
            # Pairwise tissue1_tissue2 PPINs
            diff_PPINs[f"{file.split('/')[-2]}"] = df['Protein1'].str.cat(df['Protein2'], sep='_').to_list()
        else: 
            PPINs.append([tissue_1, tissue_2, np.nan])
            PPINs.append([tissue_2, tissue_1, np.nan])
            
    
    # Pairwise tissue PPIN
    for x in list(set([x[0] for x in PPINs] + [x[1] for x in PPINs])):
        PPINs.append([x, x, 0])
    
    PPIN_distance = pd.DataFrame(PPINs, columns=['PPIN1', 'PPIN2', 'PPIN_distance']).pivot(index='PPIN1', columns='PPIN2', values='PPIN_distance')
    PPIN_similarity_mask = PPIN_distance.isnull()
    PPIN_similarity = 1 - PPIN_distance/(PPIN_distance.max().max())
    PPIN_similarity.fillna(0, inplace=True)
    
    cmap = plt.colormaps.get_cmap('viridis')
    cmap.set_bad("#f0f0f0")
    sim_plot = sns.clustermap(PPIN_similarity, mask=PPIN_similarity_mask, linewidths=0.5, figsize=(11, 10), cmap=cmap, cbar_pos=(1, .3, .03, .4))
    sim_plot.savefig(sim_plot_path, dpi=300)
    
    
    # Pairwise tissue1_tissue2 PPINs
    PPIN_similarity = [[x, y, jaccard_similarity(diff_PPINs[x], diff_PPINs[y])] for x in diff_PPINs.keys() for y in diff_PPINs.keys()]
    PPIN_similarity = pd.DataFrame(PPIN_similarity, columns=['PPIN1', 'PPIN2', 'Jaccard_similarity']).pivot(index='PPIN1', columns='PPIN2', values='Jaccard_similarity')

    sim_plot = sns.clustermap(PPIN_similarity, linewidths=0.5, figsize=(11, 10), cmap='viridis', cbar_pos=(1, .3, .03, .4))
    sim_plot.savefig(diff_sim_plot_path, dpi=300)
    

                     
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--PPICompare_path', dest='PPICompare_path', type=str, help='Path to all results PPICompare analysis.')
    parser.add_argument('--result_path', dest='result_path', type=str, help='Path to all results of PPIN_overlap analysis.')
    args = parser.parse_args()

    PPICompare_path = args.PPICompare_path
    result_path = args.result_path
    
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    get_PPIN_overlap(PPICompare_path=PPICompare_path, result_path=result_path)