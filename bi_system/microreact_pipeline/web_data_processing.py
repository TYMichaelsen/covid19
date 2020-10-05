import pandas as pd 
import logging

from utilities import datestr_to_week_func

def get_sequenced_grouped_by_week_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df[(df.genome_qc == 'MQ') | (df.genome_qc == 'HQ')]\
        .groupby(['SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'SampleDate':'week'})

def get_sequenced_grouped_by_lineage_week_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df[(df.genome_qc == 'MQ') | (df.genome_qc == 'HQ')]\
        .groupby(['SampleDate', 'nextstrain'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'SampleDate':'week', 'nextstrain':'lineage'})

def get_all_grouped_by_week_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())
    
    return df.groupby(['SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'SampleDate':'week'})

def get_sequenced_grouped_by_region_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df[(df.genome_qc == 'MQ') | (df.genome_qc == 'HQ')]\
        .groupby(['Region', 'SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'Region':'region', 'SampleDate':'week'})

def get_sequenced_grouped_by_lineage_region_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df[(df.genome_qc == 'MQ') | (df.genome_qc == 'HQ')]\
        .groupby(['Region', 'nextstrain', 'SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'Region':'region', 'nextstrain':'lineage', 'SampleDate':'week'})
    
def get_all_grouped_by_region_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df.groupby(['Region', 'SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'Region':'region', 'SampleDate':'week'})

def get_sequenced_grouped_by_age_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df[(df.genome_qc == 'MQ') | (df.genome_qc == 'HQ')]\
        .groupby(['SampleAgeGrp', 'SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'SampleDate':'week', 'SampleAgeGrp':'age_group'})

def get_sequenced_grouped_by_lineage_age_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df[(df.genome_qc == 'MQ') | (df.genome_qc == 'HQ')]\
        .groupby(['SampleAgeGrp', 'nextstrain', 'SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'nextstrain':'lineage', 'SampleDate':'week', 'SampleAgeGrp':'age_group'})
    
def get_all_grouped_by_age_df(df, config):
    df['SampleDate'] = df['SampleDate'].apply(datestr_to_week_func())

    return df.groupby(['SampleAgeGrp', 'SampleDate'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'SampleAgeGrp':'age_group', 'SampleDate':'week'})      
    
