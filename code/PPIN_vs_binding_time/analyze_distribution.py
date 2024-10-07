import pandas as pd
import numpy as np 

def roman_to_int(roman_array: list):
    roman_key = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII']
    return [roman_key.index(x) for x in roman_array]

def get_cluster_distance(annotated_PPIN: pd.DataFrame):
    print("="*20, "Get cluster distance", "="*20)
    annotated_PPIN['cluster_distance'] = abs(np.subtract(roman_to_int(annotated_PPIN['cluster_1']), roman_to_int(annotated_PPIN['cluster_2'])))
    return annotated_PPIN

def get_cluster_overlap(annotated_PPIN: pd.DataFrame):
    print("="*20, "Get cluster overlap", "="*20)
    annotated_PPIN['cluster_overlap'] = (annotated_PPIN['cluster_1'] == annotated_PPIN['cluster_2']).astype(int)
    return annotated_PPIN

def get_binding_time_difference(annotated_PPIN: pd.DataFrame):
    print("="*20, "Get binding time difference", "="*20)
    annotated_PPIN['binding_time_difference'] = abs(annotated_PPIN['Maximal_time_1'] - annotated_PPIN['Maximal_time_2'])
    return annotated_PPIN

def jaccard_similarity(list1, list2):
    list1 = list1 - set(['', 'nan'])
    list2 = list2 - set(['', 'nan'])
    if len(list1.union(list2)) == 0:
        return 0
    else :
        return(len(list1.intersection(list2)) / len(list1.union(list2)))

def get_compartment_overlap(annotated_PPIN: pd.DataFrame):
    print("="*20, "Get compartment overlap", "="*20)
    annotated_PPIN['compartment_overlap'] = annotated_PPIN.apply(lambda x: jaccard_similarity(x['Compartments_1'], x['Compartments_2']), axis=1)
    return annotated_PPIN

def plot_distribution(annotated_PPI: pd.DataFrame):
    print("="*20, "Plot distribution", "="*20)
    annotated_PPI['cluster_distance'].plot.hist(bins=6, alpha=0.5)
    return annotated_PPI

