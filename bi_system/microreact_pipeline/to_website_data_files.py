import pandas as pd
import logging
import pysftp
import json

from os import listdir
from os.path import isfile, join

from utilities import datestr_to_week_func, nut3_to_nut2_func, get_linelist
from convert_to_microreact_files import FIELD

LOGGER = logging.getLogger('to website')

def save_website_files(config, data):   
    _save_seq_grouped_by_week(data, config)
    _save_seq_grouped_by_region(data, config)
    _save_seq_grouped_by_age(data, config)

    _save_seq_grouped_by_lineage_week(data, config)
    _save_seq_grouped_by_lineage_region(data, config)
    _save_seq_grouped_by_lineage_age(data, config)

    linelist_data = get_linelist(config)
    _save_all_grouped_by_week(linelist_data, config)
    _save_all_grouped_by_region(linelist_data, config)
    _save_all_grouped_by_age(linelist_data, config)

def upload_web_files(config):
    host, username, password = config['web_host'], config['web_user'], config['web_password']
    path = config['web_path']

    file_dir = config['out_web_dir']
    files = [join(file_dir, f) for f in listdir(file_dir) if isfile(join(file_dir, f))]
    
    try:
        for file in files:
            LOGGER.debug(file)
        # srv = pysftp.Connection(host, username, password = password) 
        # with srv.cd(path):
        #     for file in files:
        #         srv.put(file)
        # srv.close()     
    except:
        LOGGER.error('Failed to sftp files to the web server.')
        # srv.close()

def _get_path(config, filename):
    return '{}/{}'.format(config['out_web_dir'], filename)

def _save_seq_grouped_by_week(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.sample_date])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={FIELD.sample_date:'week'})
     
    path = _get_path(config, 'sequenced_by_week.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _save_seq_grouped_by_lineage_week(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.sample_date, FIELD.lineage])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={FIELD.sample_date:'week', FIELD.lineage:'lineage'})

    path = _get_path(config, 'sequenced_by_lineage_week.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)
        
def _save_seq_grouped_by_region(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.region])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'region__autocolor':'region'})

    data_df = _map_region(data_df)

    path = _get_path(config, 'sequenced_by_region.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _save_seq_grouped_by_lineage_region(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.region, FIELD.lineage])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={FIELD.region:'region', FIELD.lineage:'lineage'})
    
    data_df = _map_region(data_df)

    path = _get_path(config, 'sequenced_by_lineage_region.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _save_seq_grouped_by_age(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.age_group])\
        .size()\
        .reset_index(name='cases')
    
    path = _get_path(config, 'sequenced_by_age.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _save_seq_grouped_by_lineage_age(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.age_group, FIELD.lineage])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={FIELD.lineage:'lineage'})
    
    path = _get_path(config, 'sequenced_by_lineage_age.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _save_all_grouped_by_week(linelist_data_df, config):
    linelist_data_df['Week']=linelist_data_df['SampleDate'].apply(datestr_to_week_func())
    data_df = linelist_data_df.groupby(['Week'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'Week':'week'})
    
    path = _get_path(config, 'all_by_week.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _save_all_grouped_by_region(linelist_data_df, config):
    linelist_data_df['NUTS2Code']=linelist_data_df['NUTS3Code'].apply(nut3_to_nut2_func())
    data_df = linelist_data_df.groupby(['NUTS2Code'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'NUTS2Code':'region'})

    data_df = _map_region(data_df)

    path = _get_path(config, 'all_by_region.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _save_all_grouped_by_age(linelist_data_df, config):
    data_df = linelist_data_df.groupby(['SampleAgeGrp', 'Week'])\
        .size()\
        .reset_index(name='cases')\
        .rename(columns={'SampleAgeGrp':'age_group', 'Week':'week'})      
    
    path = _get_path(config, 'all_by_age.json')
    LOGGER.info('Saving file to {}'.format(path))
    _save_df(data_df, path)

def _map_region(data_df):
    data_df['region'] = data_df['region'].map({
        'DK01':'Hovedstaden',
        'DK02':'Sjælland',
        'DK03':'Syddanmark',
        'DK04':'Midtjylland',
        'DK05':'Nordjylland'
    })
    return data_df

def _save_df(data_df, path):
    data_json = data_df.to_json(orient='records')
    data_json = json.loads(data_json)
    with open(path, 'w') as f:
        json.dump(data_json, f)
