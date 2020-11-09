import numpy as np
import pandas as pd 

def get_mink_infected_cases(df_clade_assignment, df_metadata):
    col = 'direct_aa_mutations'
    mutation = 'I692V'
    df_clade_assignment[col] = df_clade_assignment[col].replace('', np.nan)
    df_clade_assignment = df_clade_assignment.dropna(subset=[col])
    df_clade_assignment = df_clade_assignment[(df_clade_assignment[col].str.contains(mutation))]
    df = pd.merge(df_clade_assignment, df_metadata, how='inner', left_on='strain', right_on='strain')
    return df

def get_related_mink_infected_cases(df_clade_assignment, df_metadata):
    col = 'direct_aa_mutations'
    df_clade_assignment[col] = df_clade_assignment[col].replace('', np.nan)
    df_clade_assignment = df_clade_assignment.dropna(subset=[col])
    df_clade_assignment = df_clade_assignment[(df_clade_assignment[col].str.contains('Y453F') & df_clade_assignment[col].str.contains('D614G'))]
    df = pd.merge(df_clade_assignment, df_metadata, how='inner', left_on='strain', right_on='strain')
    return df
