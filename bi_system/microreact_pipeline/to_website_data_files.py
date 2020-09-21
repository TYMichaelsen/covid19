import pandas as pd
import logging

from utilities import datestr_to_week_func, nut3_to_nut2_func, get_linelist
from convert_to_microreact_files import FIELD

LOGGER = logging.getLogger('to website')

def save_website_files(config, data):   
    _save_seq_grouped_by_week(data, config)
    _save_seq_grouped_by_region(data, config)
    _save_seq_grouped_by_age(data, config)

    linelist_data = get_linelist(config)
    _save_all_grouped_by_week(linelist_data, config)
    _save_all_grouped_by_region(linelist_data, config)
    _save_all_grouped_by_age(linelist_data, config)

def _get_path(config, filename):
    return '{}/{}'.format(config['out_web_dir'], filename)

def _save_seq_grouped_by_week(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.epi_week]).size().reset_index(name='cases')
    
    path = _get_path(config, 'sequenced_by_week.csv')
    LOGGER.info('Saving file to {}'.format(path))
    data_df.to_csv(path)
        
def _save_seq_grouped_by_region(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.region]).size().reset_index(name='cases')
    
    path = _get_path(config, 'sequenced_by_region.csv')
    LOGGER.info('Saving file to {}'.format(path))
    data_df.to_csv(path)

def _save_seq_grouped_by_age(data, config):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([FIELD.age_group]).size().reset_index(name='cases')
    
    path = _get_path(config, 'sequenced_by_age.csv')
    LOGGER.info('Saving file to {}'.format(path))
    data_df.to_csv(path)

def _save_all_grouped_by_week(linelist_data_df, config):
    linelist_data_df['Week']=linelist_data_df['SampleDate'].apply(datestr_to_week_func())
    data_df = linelist_data_df.groupby(['Week']).size().reset_index(name='cases')
    
    path = _get_path(config, 'all_by_week.csv')
    LOGGER.info('Saving file to {}'.format(path))
    data_df.to_csv(path)

def _save_all_grouped_by_region(linelist_data_df, config):
    linelist_data_df['NUTS2Code']=linelist_data_df['NUTS3Code'].apply(nut3_to_nut2_func())
    data_df = linelist_data_df.groupby(['NUTS2Code']).size().reset_index(name='cases')

    path = _get_path(config, 'all_by_region.csv')
    LOGGER.info('Saving file to {}'.format(path))
    data_df.to_csv(path)

def _save_all_grouped_by_age(linelist_data_df, config):
    data_df = linelist_data_df.groupby(['SampleAgeGrp']).size().reset_index(name='cases')       
    
    path = _get_path(config, 'all_by_age.csv')
    LOGGER.info('Saving file to {}'.format(path))
    data_df.to_csv(path)