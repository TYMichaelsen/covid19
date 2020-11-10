import sys
import json
import pysftp
import os

from datetime import datetime
from pandas import read_csv, option_context

def get_linelist(config):

    path = ''
    for f in os.scandir(config['raw_ssi_file']):
        if 'linelist' in f.name:
            path = '{}/{}'.format(config['raw_ssi_file'], f.name)
            break
    
    assert path != '', 'linelist file not found.'
    linelist = read_csv(path, sep='\t')
    assert linelist.empty == False, 'failed to load linelist file'   
    return linelist

def _get_latest_linelist_file(config):
    path = '/'.join(config['raw_ssi_file'].split('/')[:-1])
    file = config['raw_ssi_file'].split('/')[-1]
    filename = file.split('_')[1]
    date_format = file.split('_')[0]
    
    dates = []
    for f in os.scandir(path):
        if filename not in f.name:
            continue
        date_str = f.name.split('_')[0]
        dates.append(datetime.strptime(date_str, date_format))
    
    assert len(dates) > 0, 'failed to find stat file.'
    latest_date_str = max(dates).strftime(date_format)
    return '{}/{}_{}'.format(path, latest_date_str, filename)

def datestr_to_week_func():
      return lambda date: datetime.strptime(date, '%Y-%m-%d').isocalendar()[1]

def datestr_to_year_func():
      return lambda date: datetime.strptime(date, '%Y-%m-%d').year

def datestr_to_week_and_year_func():
      return lambda date: '{},{}'.format(
            datetime.strptime(date, '%Y-%m-%d').isocalendar()[1],
            datetime.strptime(date, '%Y-%m-%d').year)

def nut3_to_nut2_func():
      return  lambda nut3: nut3[:-1] if type(nut3) is str else ''

def save_df_as_json(data_df, path):
    data_json = data_df.to_json(orient='records')
    data_json = json.loads(data_json)
    with open(path, 'w') as f:
        json.dump(data_json, f)

def stfp_file(host, user, password, destination_path, local_path, logger):
    try:
        logger.info('Attempting to connect to {}'.format(host))
        srv = pysftp.Connection(host, user, password=password)
        logger.info('Connected to {}'.format(host))
        with srv.cd(destination_path):
            logger.info('Moved to {}'.format(destination_path))
            logger.info('Sending {}...'.format(local_path))
            srv.put(local_path, preserve_mtime=True)
        srv.close()
    except:
        e = sys.exc_info()
        logger.error('Failed to sftp {} to the web server.'.format(local_path))
        logger.error(e)
        srv.close()

def print_full(df):
      with option_context('display.max_rows', None, 'display.max_columns', None):
            print(df)
