import pandas as pd
from itertools import chain


def clean_up_compartment_terms(compartments: str, sep=','):
    return list(compartments.replace(' ', '').replace('noprediction', '').replace('N/A', '').replace('-', '').replace('miscellaneous', '').lower().split(sep))

def get_compartment_list(mRNPchrono_table: pd.DataFrame):
    print("="*20, "Get compartment list", "="*20)
    compartment_list = mRNPchrono_table.apply(lambda x: [clean_up_compartment_terms(str(x['local_HPA']), sep=';'), 
                                                         clean_up_compartment_terms(str(x['local_HCM_NMF'])),
                                                         clean_up_compartment_terms(str(x['local_HCM_SAFE']))], 
                                              axis=1)
    mRNPchrono_table['Compartments'] = [set(chain.from_iterable(x)) for x in compartment_list]
    return(mRNPchrono_table)


def annotate_by_mRNPchrono(PPIN: pd.DataFrame, mRNPchrono_table: pd.DataFrame, cleanup: bool = True):
    print("="*20, "Annotate PPIN", "="*20)
    print(PPIN.head())
    print(mRNPchrono_table.head())
    annotated_PPIN = PPIN.merge(mRNPchrono_table, left_on='Protein1', right_on='MasterAccession', how='left')
    annotated_PPIN = annotated_PPIN.merge(mRNPchrono_table, left_on='Protein2', right_on='MasterAccession', how='left',
                                         suffixes=('_1', '_2'))
    annotated_PPIN.columns = annotated_PPIN.columns.str.replace(' ', '_')
    
    if (cleanup):
        annotated_PPIN = clean_annotated_PPIN(annotated_PPIN)
    return annotated_PPIN


def clean_annotated_PPIN(annotated_PPIN: pd.DataFrame):
    print("="*20, "Clean annotated PPIN", "="*20)
    print("Before cleaning: ", annotated_PPIN.shape)
    cleaned_PPIN = annotated_PPIN.dropna(subset=['MasterAccession_1', 'MasterAccession_2', 'cluster_1', 'cluster_2'])
    print("After cleaning: ", cleaned_PPIN.shape)
    return cleaned_PPIN