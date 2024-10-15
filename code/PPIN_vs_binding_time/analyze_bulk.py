#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import starbars
from scipy import stats

if __name__ == '__main__':  
    sns.set_style("whitegrid")
    parser = argparse.ArgumentParser()
    parser.add_argument('--PPIN_path', dest='PPIN_path', type=str, help='Path to all results of PPIN_vs_binding_time analysis.')
    parser.add_argument('--plot_path', dest='plot_path', type=str, help='Path to  of PPIN_vs_binding_time analysis.')
    args = parser.parse_args()

    PPIN_path = args.PPIN_path
    plot_path = args.plot_path
    tissue = PPIN_path.split('/')[-1]

    path_to_sample_PPIN = [file for file in os.listdir(PPIN_path) if file.startswith("annotated_")]
    tissue_specific_PPIN = [pd.read_csv(os.path.join(PPIN_path, x), sep='\t') for x in path_to_sample_PPIN]
    tissue_specific_PPIN = pd.concat(tissue_specific_PPIN)

    
    fig = plt.figure(figsize=(8, 6)) 
    
    ax3 = fig.add_subplot(2, 2, 1)
    sns.boxplot(data=tissue_specific_PPIN, x="sample", y='binding_time_difference', hue="is_PPI", hue_order=[0, 1],
                palette='plasma', gap=0.1, legend=False)
    ax3.set_title('Binding time difference')
    ax3.set_xlabel('Samples')
    ax3.set_ylabel('Binding time difference (min)')
    ax3.set_ylim(0, 250)
    # legend = ax3.get_legend()
    # handles = legend.legend_handles
    # ax3.legend(handles, ['Non-PPIN', 'PPIN'], title='')
    # annotations = [(('1', '0'), ('2', '0')), (('1', '1'), ('2', '1'))]
    # for sample in [1, 2]:
    #     a = tissue_specific_PPIN.query(f'is_PPI == 0 & sample == {sample}')['binding_time_difference']
    #     b = tissue_specific_PPIN.query(f'is_PPI == 1 & sample == {sample}')['binding_time_difference']
    #     annotations[sample - 1] = annotations[sample - 1] + (stats.mannwhitneyu(a, b)[1], )
    # starbars.draw_annotation(annotations)
    
    
    ax2 = fig.add_subplot(2, 2, 2)
    for sample in tissue_specific_PPIN['sample'].unique():
        sns.histplot(data=tissue_specific_PPIN[tissue_specific_PPIN['sample']==sample], x="binding_time_difference", hue="is_PPI", hue_order=[0, 1], 
                     palette='plasma', bins=250, multiple="layer", shrink=0.8, stat='percent', common_norm=False, alpha=0.5, 
                     element='step', fill=False, cumulative=False, legend=True)
    ax2.set_title(f'Binding time difference\ncumulative distribution')
    ax2.set_xlabel('Binding time difference')
    legend = ax2.get_legend()
    handles = legend.legend_handles
    ax2.legend(handles, ['Non-PPIN', 'PPIN'], title='')

    
    ax1 = fig.add_subplot(2, 2, 3)
    for sample in tissue_specific_PPIN['sample'].unique():
        sns.histplot(data=tissue_specific_PPIN[tissue_specific_PPIN['sample']==sample], x="cluster_distance", hue="is_PPI", hue_order=[0, 1],
                     palette='plasma', bins=7, multiple="layer", shrink=0.8, stat='percent', common_norm=False, alpha=0.5, 
                     element='step', fill=False, cumulative=True, legend=True)
    ax1.set_title('Cluster distance distribution')
    ax1.set_xlabel('Cluster distance')
    ax1.set_ylim(0, 100)
    legend = ax1.get_legend()
    handles = legend.legend_handles
    ax1.legend(handles, ['Non-PPIN', 'PPIN'], title='')
    
    
    ax4 = fig.add_subplot(2,2,4)
    for sample in tissue_specific_PPIN['sample'].unique():
        sns.histplot(data=tissue_specific_PPIN[tissue_specific_PPIN['sample']==sample], x="cluster_overlap", hue="is_PPI", palette='plasma', 
                    hue_order=[0, 1], bins=2, multiple="dodge", shrink=0.995, stat='percent', 
                    common_norm=False, alpha=0.5, fill=True, cumulative=False, legend=True)
    ax4.set_title('Cluster overlap')
    ax4.set_xlabel('Cluster overlap')
    ax4.set_xticks([True, False])
    ax4.set_xticklabels(['Overlap', 'No overlap']) 
    ax4.set_ylim(0, 100)
    legend = ax4.get_legend()
    handles = legend.legend_handles
    ax4.legend(handles, ['Non-PPIN', 'PPIN'], title='')
    
    fig.tight_layout()
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    