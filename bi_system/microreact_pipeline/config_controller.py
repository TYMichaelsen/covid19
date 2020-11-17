import os
import json
import pathlib
import logging
import re

from datetime import datetime
from data_cleansing_metadata import check_file
from shutil import copyfile

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

def set_config_paths(config):
    
    config = set_linelist_path(config)
    paths = get_nextstrain_paths(config)
    for base_path in paths:
        try:
            m_path = get_metadata_path(base_path, config)
            c_path = get_clade_path(base_path, config)
            t_path = get_tree_path(base_path, config)
            
            for fileame in [m_path, c_path, t_path]:
                if os.path.isfile(fileame) == False:
                    raise Exception('{} is not a file. Retrying earlier date.')

            config['sequenced_metadata_file'] = m_path
            config['global_clade_assignment_file'] = c_path
            config['clade_tree_path'] = t_path
            return config
        except:
            continue
    raise Exception('Did not find consistent date for metadata, tree and clade assignment.')

def set_linelist_path(config):
    base_path = config['raw_ssi_file']
    with os.scandir(base_path) as it:
        for file in it:
            if 'linelist' in file.name and file.name[-1] != '#':
                config['raw_ssi_file'] = '{}/{}'.format(base_path, file.name)
                return config
    raise Exception('Unable to find linelist file.')

def get_tree_path(base_path, config):
    config_path = config['clade_tree_path']
    config_path = config_path.split('/nextstrain/')[1]
    config_path = config_path.split('/')[1:]
    config_path = '/'.join(config_path)
    path = '{}/{}'.format(base_path,config_path)
    return path

def get_clade_path(base_path, config):
    config_path = config['global_clade_assignment_file']
    config_path = config_path.split('/nextstrain/')[1]
    config_path = config_path.split('/')[1:]
    config_path = '/'.join(config_path)
    path = '{}/{}'.format(base_path,config_path)
    return path

def get_metadata_path(base_path, config):
    config_path = config['sequenced_metadata_file']
    config_path = config_path.split('/nextstrain/')[1]
 
    file_name = config_path.split('/')[-1]
    config_path = config_path.split('/')[1:-1]
    config_path = '/'.join(config_path)
    file_suffix = '_'.join(file_name.split('_')[1:])
    path = '{}/{}'.format(base_path,config_path)
    
    for file in os.scandir(path):
        if file_suffix not in file.name:
            continue
        return '{}/{}'.format(path, file.name)
    raise Exception('Unable to find metadata file.')

def get_nextstrain_paths(config):
    path = get_nextstrain_base_path(config)
    path_prefix = path.split('/nextstrain/')[0] + '/nextstrain'
    path_suffix = path.split('/nextstrain/')[1]
    date_format = path_suffix.split('/')[0].split('_')[0]
    dates = get_nextstrain_date_list(path_prefix, date_format)
    paths = get_paths_from_date_list(dates, date_format, path_prefix, path_suffix)
    return paths
        
def get_paths_from_date_list(dates, date_format, path_prefix, path_suffix):
    paths = []
    for date in dates:
        dir_prefix = date.strftime(date_format)
        dir_suffix = path_suffix.split('/')[0].split('_')[1]
        dir_name = '{}_{}'.format(dir_prefix, dir_suffix)
        path = '{}/{}'.format(path_prefix, dir_name)
        paths.append(path)
    return paths

def get_nextstrain_base_path(config):
    for k in config.keys():
        if 'nextstrain' not in config[k]: 
            continue
        return config[k]
    raise Exception('Unable to find nextstrain path in config.')

def get_nextstrain_date_list(path, date_format):
    dates = []
    for date_dir in os.scandir(path):
        if date_dir.is_dir() == False:
            continue
        if date_dir.name == 'latest':
            continue                       
        try:
            date_str = date_dir.name.split('_')[0]
            date = datetime.strptime(date_str, date_format)
            dates.append(date)
        except:
            LOGGER.warning('Unrecognized date format of date string: {}. Continuing...'.format(date_str))
    return sorted(dates, reverse=True)