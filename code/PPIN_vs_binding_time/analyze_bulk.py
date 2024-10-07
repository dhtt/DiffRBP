import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from os import listdir, path
import starbars
from scipy import stats

if __name__ == '__main__':
    tissue = 'adipose'
    path_to_annotated_PPINs = f"../../data/RBP2GO/PPIN_vs_binding_time/{tissue}/"
    path_to_sample_PPIN = [file for file in listdir(path_to_annotated_PPINs) if file.startswith("annotated_")]
    tissue_specific_PPIN = [pd.read_csv(path.join(path_to_annotated_PPINs, x), sep='\t') for x in path_to_sample_PPIN]
    tissue_specific_PPIN = pd.concat(tissue_specific_PPIN)
    print(tissue_specific_PPIN.head())
    print(tissue_specific_PPIN.tail())
    
    
    g = sns.FacetGrid(tissue_specific_PPIN, col="is_PPI", margin_titles=True, hue="sample", palette='flare', height=4)
    g.map(sns.histplot, "cluster_distance", bins=7, multiple="dodge", shrink=0.8, stat='percent', 
          alpha=0.7, element='step', fill=False, cumulative=False, legend=True)
    g.add_legend(title="Sample")
    # axes = g.axes.flatten()
    # axes[1].set_title("PPIN")
    # axes[0].set_title("Non-PPIN")
    g.set(ylim=(0, 100), xlabel='Cluster distance')
    g.savefig(path.join(path_to_annotated_PPINs, f"cluster_distance.png"))
    
    
    g = sns.FacetGrid(tissue_specific_PPIN, col="is_PPI", margin_titles=True, hue="sample", palette='flare', height=4, sharey=False)
    g.map(sns.histplot, "cluster_overlap", multiple="dodge", stat='count', alpha=0.7, legend=True)
    g.add_legend(title="Sample")
    g.set(xlabel='Cluster overlap', xticks=[0, 1], xticklabels=['No', 'Yes'])
    g.savefig(path.join(path_to_annotated_PPINs, f"cluster_overlap.png"))
    
    
    g = sns.FacetGrid(tissue_specific_PPIN, col="is_PPI", margin_titles=True, hue="sample", palette='flare', height=4)
    g.map(sns.histplot, "binding_time_difference", bins=15, multiple="dodge", shrink=0.8, stat='percent', 
          alpha=0.7, element='step', fill=False, cumulative=False, legend=True)
    g.add_legend(title="Sample")
    g.set(ylim=(0, 100), xlabel='Binding Time Difference')
    g.savefig(path.join(path_to_annotated_PPINs, f"binding_time_difference_cumulative.png"))
    
    
    plt1, ax1 = plt.subplots()
    ax1 = sns.boxplot(data=tissue_specific_PPIN, hue='sample', x='is_PPI', y='binding_time_difference', 
            palette='flare', gap=0.1)
    
    annotations = [(('1', '0'), ('1', '1')), (('2', '0'), ('2', '1'))]
    for sample in [1, 2]:
        a = tissue_specific_PPIN.query(f'is_PPI == 0 & sample == {sample}')['binding_time_difference']
        b = tissue_specific_PPIN.query(f'is_PPI == 1 & sample == {sample}')['binding_time_difference']
        annotations[sample - 1] = annotations[sample - 1] + (stats.mannwhitneyu(a, b)[1], )
    starbars.draw_annotation(annotations)
    
    plt1.savefig(path.join(path_to_annotated_PPINs, f"binding_time_difference_box.png"))