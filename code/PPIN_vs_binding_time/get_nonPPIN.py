import pandas as pd
from annotate_by_mRNPchrono import *
from itertools import combinations

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