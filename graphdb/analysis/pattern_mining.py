import pandas as pd
import math 
import numpy as np 
import matplotlib.pyplot as plt
from datetime import datetime

def get_df(path):
    df = pd.read_csv('temps/linelist.tsv', sep='\t')
    df_municipalities = pd.read_csv('../../bi_system/stable_dims/municipalities.tsv', sep='\t')
    df_clades = pd.read_csv('temps/clade_assignment.tsv', sep='\t')
    
    # INNER JOIN HERE REMOVES 54K / 66K ENTRIES
    df = pd.merge(df, df_municipalities, how='inner', left_on='MunicipalityCode', right_on='code')
    df = pd.merge(df, df_clades, how='inner', left_on='ssi_id', right_on='strain')

    df['date_linelist'] = df['date_linelist'].apply(lambda x:datetime.strptime(x, '%Y-%m-%d'))
    df = df[['strain','code', 'name', 'administrative_center', 'area', 'population', 'Region','zipcode_name','NUTS3Text','Address', 'clade_y', 'cluster_core', 'cluster_size', 'cluster_date', 'direct_mutations', 'direct_aa_mutations', 'date_linelist']]
    df = df.rename(columns={
        'strain':'id',
        'code':'municipality_code',
        'name':'municipality_name',
        'Region':'region',
        'NUTS3Text':'nuts3_name',
        'Address':'address',
        'clade_y':'clade'
    })    
    return df

def filter_df(df, cap):
    df_grp = df.groupby(['municipality_name']).size().reset_index(name='cases')
    municipality_blacklist = df_grp[df_grp['cases'] < cap]['municipality_name']
    df = df[~df['municipality_name'].isin(municipality_blacklist)]
    return df

def get_top_value(lst):
    vals, freq = np.unique(lst, return_counts=True)
    argmax = np.argmax(freq)
    return vals[argmax]

def create_dims(df):
    df['direct_mutations'] = df['direct_mutations'].apply(lambda x: np.array(x.split(',')) if str(x) != 'nan' else [])   
    df['direct_aa_mutations'] = df['direct_aa_mutations'].apply(lambda x: np.array(x.split(';')) if (str(x) != 'nan') else [])

    df['no_direct_mutations'] = df['direct_mutations'].apply(lambda x: len(x))
    df['no_direct_aa_mutations'] = df['direct_aa_mutations'].apply(lambda x: len(x))
    df['clade_depth'] = df['clade'].apply(lambda x: len(x.split('/')) if str(x) != 'nan' else 0)
    
    #AGGR
    df_agg = df.groupby(['municipality_name']).agg({
        'clade_depth':['mean'],
        'no_direct_mutations':['mean'],

    })
    df_agg.columns = ["_".join(x) for x in df_agg.columns.ravel()]

    #CASES
    df_cases = df.groupby(['municipality_name']).size().reset_index(name='no_cases')
   
    #CLADES
    df_clades = df.groupby(['municipality_name'])['clade'].apply(list).reset_index(name='clades')
    df_clades['no_unique_clades'] = df_clades['clades'].apply(lambda x: len(np.unique(np.array(x))))
    df_clades['top_clade'] = df_clades['clades'].apply(get_top_value)
    df_clades = df_clades[['municipality_name', 'no_unique_clades', 'top_clade']]

    #MUTATIONS
    df_mutations = df.groupby(['municipality_name'])['direct_mutations'].apply(list).reset_index(name='mutations')
    df_mutations['mutations'] = df_mutations['mutations'].apply(lambda x: np.concatenate(np.array(x)))
    df_mutations['no_unique_direct_mutations'] = df_mutations['mutations'].apply(lambda x: len(np.unique(x)))
    df_mutations['top_direct_mutation'] = df_mutations['mutations'].apply(get_top_value)
    df_mutations = df_mutations[['municipality_name', 'no_unique_direct_mutations', 'top_direct_mutation']]

    #AA MUTATIONS
    df_aa_mutations = df.groupby(['municipality_name'])['direct_aa_mutations'].apply(list).reset_index(name='mutations')
    df_aa_mutations['mutations'] = df_aa_mutations['mutations'].apply(lambda x: np.concatenate(np.array(x)))
    df_aa_mutations['no_unique_direct_aa_mutations'] = df_aa_mutations['mutations'].apply(lambda x: len(np.unique(x)))
    df_aa_mutations['top_direct_aa_mutation'] = df_aa_mutations['mutations'].apply(get_top_value)
    df_aa_mutations = df_aa_mutations[['municipality_name', 'no_unique_direct_aa_mutations', 'top_direct_aa_mutation']]

    #PATIENTS FROM THE BEGINNING OF SECOND WAVE (2020-08-02)
    print(len(df.index))
    df_patients = df[df['date_linelist'] > datetime.strptime('2020-08-02', '%Y-%m-%d')]
    print(len(df_patients))
    # df_patients = df.groupby('date_linelist').size().reset_index(name='cases')
    

    # df_patients.plot.bar(x='date_linelist', y='cases')    
    


    # df = pd.merge(df_agg, df_cases, how='inner', on='municipality_name')
    # df = pd.merge(df, df_clades, how='inner', on='municipality_name')
    # df = pd.merge(df, df_mutations, how='inner', on='municipality_name')
    # df = pd.merge(df, df_aa_mutations, how='inner', on='municipality_name')

    # print(df)
    
    

if __name__ == '__main__':
    df = get_df('temps/linelist.tsv')
    # df = filter_df(df, 10)


    create_dims(df)
    
