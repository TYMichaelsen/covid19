import logging
import pysftp
import sys

from pandas import read_csv

LOGGER = logging.getLogger('website stat files')

def send_stats(config):
    received, sequenced, genomes, uploaded = _load_stats_from_file(
        config['stats_file'])
    stat_dict = {
        'received': received,
        'sequenced': sequenced,
        'genomes': genomes,
        'uploaded': uploaded
    }
    # _stfp_stats(config, stat_dict)

def _stfp_stats(config, data):
    try:
        srv = pysftp.Connection(config['host'], config['user'], config['password'])
        with srv.cd(config['web_stats_path']):
            LOGGER.debug('Sending stats...')
            srv.put(data)
        srv.close()
    except:
        e = sys.exc_info()
        LOGGER.error('Failed to sftp stats to the web server.')
        LOGGER.error(e)
        srv.close()

def _load_stats_from_file(path):
    df = read_csv(path, sep="\t")
    return df['Recieved'][0], df['Sequenced'][0], df['Genomes'][0], df['Uploaded'][0]
