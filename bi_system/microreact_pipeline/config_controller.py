import os
import json
import pathlib
import logging
import re

from datetime import datetime
from data_cleansing_metadata import check_file

LOGGER = logging.getLogger("config controller")

def get_config(config_filepath, config_template_filepath):
    config_filepath = check_file(config_filepath)
    with open(config_filepath, 'r') as f:
        config = json.load(f)
    
    with open(config_template_filepath, 'r') as f:
        expected_config = json.load(f)
    
    for k in expected_config.keys():
        if k not in config.keys():
            err = 'Error: expected {} in config file {} but could not find such a key.'.format(k, config_filepath)
            LOGGER.error(err)
            raise KeyError(err)
    return config

def set_config_nextstrain(config, set_date_str='latest', date_folder_suffix='nextstrain'): 
    modified_config = {}
    for k in config.keys():
        if 'nextstrain' not in config[k]: 
            modified_config[k] = config[k]
            continue

        replacement_path = _get_replacement_nextstrain_path(config[k], set_date_str, date_folder_suffix)
        modified_config[k] = replacement_path
    return modified_config
   
def set_config_linelist(config, option):
    path = config['raw_ssi_file']
    if option != 'latest':
        LOGGER.info('Linelist is set to {}'.format(path))
        return config
    
    replacement_path = _get_replacement_linelist_path(path)
    LOGGER.info('Linelist is set to {}'.format(replacement_path))
    config['raw_ssi_file'] = replacement_path
    return config
    
def update_latest_nextstrain(config):
    common_path_prefix = _get_next_strain_path(config)
    nextstrain_dates = sorted(_get_nextstrain_date_list(common_path_prefix), reverse=True)

    for date in nextstrain_dates:
        try:
            date_str = _get_datestr_from_date(date)
            for k in config.keys():
                if 'nextstrain' not in config[k]:
                    continue                
                
                _, _, suffix = _split_path(config[k])
                source_path = '{}/{}_nextstrain/{}'.format(common_path_prefix, date_str, suffix)
                target_path = '{}/latest/{}'.format(common_path_prefix, suffix)
                
                LOGGER.info("Copying {} to latest...".format(source_path))
                if os.path.exists(os.path.dirname(target_path)) == False:
                    os.makedirs(os.path.dirname(target_path))
                copyfile(source_path, target_path)
            LOGGER.info("Done.")
            return
        except:
            LOGGER.warning("Unable to copy {} to latest. Retrying earlier date...".format(source_path))

def _get_replacement_nextstrain_path(path, date_str, date_folder_suffix):
    prefix, _, suffix = _split_path(path)
    if date_str == 'latest':
        return '{}/latest/{}'.format(prefix, suffix)

    if os.path.isdir('{}/{}_{}'.format(prefix, date_str, date_folder_suffix)):
        return '{}/{}_{}/{}'.format(prefix, date_str, date_folder_suffix, suffix)

    dir_str = '{}/{}_{}'.format(prefix, date_str, date_folder_suffix)
    raise NotADirectoryError('Directory {} does not exist, make sure to provide correct date and date suffix.'.format(dir_str))

def _get_nextstrain_date_list(path):
    dates = []
    for date_dir in os.scandir(path):
        if date_dir.is_dir() == False:
            continue
        if date_dir.name == 'latest':
            continue               
        date_str = date_dir.name.split('_')[0]
        date_format='%Y-%m-%d-%H-%M' if len(date_str.split('-')) > 3 else '%Y-%m-%d'

        try:
            date = datetime.strptime(date_str, date_format)
            dates.append(date)
        except:
            LOGGER.error('Unrecognized date format of date string: {}'.format(date_str))
    return dates

def _split_path(path):
    path_prefix = path.split('/nextstrain/')[0] + '/nextstrain'
    path_suffix = pathlib.Path(path.split('/nextstrain/')[1])
    path_suffix_after_date = str(pathlib.Path(*path_suffix.parts[1:]))
    date_dir = str(pathlib.Path(path_suffix.parts[0]))
    return path_prefix, date_dir, path_suffix_after_date

def _get_next_strain_path(config):
    for k in config.keys():   
        if 'nextstrain' in config[k]:
            return _split_path(config[k])[0]
    raise Exception('Unable to identify nextstrain path in the config file.')

def _get_datestr_from_date(date):
    if date.hour == 0 and date.minute == 0:
        return date.strftime('%Y-%m-%d')
    return date.strftime('%Y-%m-%d-%H-%M')

def _get_linelist_date_list(path):
    dates = []
    pattern = re.compile("^Lineliste_\\d{6}.xlsx$")
    
    for file in os.scandir(path):
        if file.is_dir():
            continue
        if pattern.match(file.name) == None:
            continue
        
        date_str = file.name.replace('_','.').split('.')[1]
        date = datetime.strptime(date_str, '%d%m%y')
        dates.append(date)
    return dates

def _get_replacement_linelist_path(path):
    linelist_dir = path.split('raw-ssi-metadata/')[0] + 'raw-ssi-metadata/'
    dates = _get_linelist_date_list(linelist_dir)
    
    max_date = max(dates)
    max_date_str = max_date.strftime('%d%m%y')
    max_filename = 'Lineliste_{}.xlsx'.format(max_date_str)
    max_path = '{}{}'.format(linelist_dir, max_filename)
    return max_path