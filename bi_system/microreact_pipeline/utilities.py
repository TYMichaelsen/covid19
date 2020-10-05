import sys
import json
import pysftp

from datetime import datetime
from pandas import read_excel

def get_linelist(config):
      linelist = read_excel(config['raw_ssi_file'])
      assert linelist.empty == False
      return linelist

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
            srv.put(local_path)
        srv.close()
    except:
        e = sys.exc_info()
        logger.error('Failed to sftp {} to the web server.'.format(local_path))
        logger.error(e)
        srv.close()