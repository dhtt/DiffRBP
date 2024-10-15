#!/usr/bin/env python3

import pandas as pd
from annotate_by_mRNPchrono import *
from itertools import combinations
import pickle
import os

def check_interaction(interaction, all_possible_interactions: dict):
    # Swapping is crucially not necessary, as the order of proteins in the interaction is already sorted
    # protein_1 = interaction['Protein1']
    # protein_2 = interaction['Protein2']
    # if protein_1 > protein_2:
    #     protein_1_, protein_2, protein_1 = protein_2, protein_1, protein_1_
        
    if not (interaction['Protein1'] in all_possible_interactions.keys()):
        return False
    else:
        if interaction['Protein2'] in all_possible_interactions[interaction['Protein1']]:
            return True
        else:
            return False


def filter_impossible_interactions(all_interactions_path: str, PPIXpress_path: str, nonPPIN_df: pd.DataFrame):
    print("="*20, "Filter impossible interactions", "="*20)
    if not os.path.exists(all_interactions_path):
        PPIN_files = []
        for root, dirs, files in os.walk(PPIXpress_path):
            for file in files:
                if file.endswith('_ppin.txt') & (not file.startswith('.')):
                    PPIN_files.append(os.path.join(root, file))
        
        PPIN = pd.concat([pd.read_csv(x, sep=' ', header=0).iloc[:, :2] for x in PPIN_files])
        PPIN = PPIN.groupby('Protein1', as_index=False).agg(lambda x: set(x))
        all_possible_interactions = dict(zip(PPIN['Protein1'], PPIN['Protein2']))
            
        with open(all_interactions_path, 'wb') as f:
            pickle.dump(all_possible_interactions, f, pickle.HIGHEST_PROTOCOL)
            
        filter_impossible_interactions(all_interactions_path, PPIXpress_path, nonPPIN_df)
        
    else:
        with open(all_interactions_path, 'rb') as f:
            all_possible_interactions = pickle.load(f)
        
        # print("BEFORE FILTERING: ", nonPPIN_df.shape)
        nonPPIN_df['is_possible'] = nonPPIN_df.apply(lambda x: check_interaction(x, all_possible_interactions), axis=1)
        # print("AFTER FILTERING: ", nonPPIN_df.shape)
        
    return nonPPIN_df


def get_nonPPIN(PPIN: pd.DataFrame, mRNPchrono_table: pd.DataFrame, cleanup: bool = True):
    print("="*20, "Get non-PPIN", "="*20)
    source_set = set(PPIN['Protein1'])
    target_set = set(PPIN['Protein2'])
    all_interactions = set(combinations(sorted(source_set.union(target_set)), 2))
    print("Number of possible interactions (annotated proteins only): ", len(all_interactions))
    
    
    interaction_set = set()
    PPIN.apply(lambda x: 
        interaction_set.add((x['Protein1'], x['Protein2'])) if x['Protein1'] < x['Protein2'] 
        else interaction_set.add((x['Protein2'], x['Protein1'])), 
        axis=1)
    print("Number of interactions (annotated proteins only): ", len(interaction_set))
    print("Overlap rate with possible interactions: ", 
          100*len(interaction_set)/len(interaction_set.intersection(all_interactions)))
    
    
    non_interaction_set = all_interactions - interaction_set
    print("Number of non-interactions (annotated proteins only): ", len(non_interaction_set))   
    nonPPIN = pd.DataFrame(list(non_interaction_set), columns=['Protein1', 'Protein2'])
    nonPPIN['weight'] = 0
    nonPPIN = annotate_by_mRNPchrono(PPIN=nonPPIN, mRNPchrono_table=mRNPchrono_table, cleanup=cleanup)
    return nonPPIN