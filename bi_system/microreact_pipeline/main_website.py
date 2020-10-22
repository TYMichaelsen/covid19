import argparse
import logging
import os

from datetime import datetime
from pandas import read_csv
from to_website_data_files import save_website_files, upload_web_files
from config_controller import get_config
from utilities import save_df_as_json, stfp_file, datestr_to_week_and_year_func
from web_data_processing import get_all_grouped_by_age_df, get_sequenced_grouped_by_lineage_age_df, get_sequenced_grouped_by_age_df, get_all_grouped_by_region_df, get_sequenced_grouped_by_lineage_region_df, get_sequenced_grouped_by_region_df, get_sequenced_grouped_by_week_df, get_sequenced_grouped_by_lineage_week_df, get_all_grouped_by_week_df

def set_logging(config):
    logging.basicConfig(level=logging.DEBUG, filename=config['log_filepath'], filemode='w')
    formatter = logging.Formatter('%(name)-12s: %(levelname)-6s %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


def _aggregate_data(config):
    logger = logging.getLogger("web data")
    
    logger.info('Opening data file: {}'.format(config['data_file']))
    df = read_csv(config['data_file'], sep="\t")
    data = []

    data.append({'name':'seq_groupedby_week', 'df':get_sequenced_grouped_by_week_df(df.copy(), config)})
    data.append({'name':'seq_groupedby_age', 'df':get_sequenced_grouped_by_age_df(df.copy(), config)})
    data.append({'name':'seq_groupedby_region', 'df':get_sequenced_grouped_by_region_df(df.copy(), config)})

    data.append({'name':'all_groupedby_week', 'df':get_all_grouped_by_week_df(df.copy(), config)})
    data.append({'name':'all_groupedby_age', 'df':get_all_grouped_by_age_df(df.copy(), config)})
    data.append({'name':'all_groupedby_region', 'df':get_all_grouped_by_region_df(df.copy(), config)})
    
    data.append({'name':'seq_groupedby_lineage_week', 'df':get_sequenced_grouped_by_lineage_week_df(df.copy(), config)})
    data.append({'name':'seq_groupedby_lineage_age', 'df':get_sequenced_grouped_by_lineage_age_df(df.copy(), config)})
    data.append({'name':'seq_groupedby_lineage_region', 'df':get_sequenced_grouped_by_lineage_region_df(df.copy(), config)})
    return data
      
def _send_stats(config):
    logger = logging.getLogger('web stats')
    filepath = _get_latest_stats_file(config)
    logger.info('Found latest stat file at: {}'.format(filepath))
    df = read_csv(filepath, sep="\t")
    save_df_as_json(df, config['stats_save_path'])
    stfp_file(config['host'], config['user'], config['password'], config['stats_server_path'], config['stats_save_path'], logger)

def _send_data(config, data):
    logger = logging.getLogger('web data')
    for e in data:
        path = '{}/{}.json'.format(config['data_save_path'], e['name'])
        save_df_as_json(e['df'], path)
        stfp_file(config['host'], config['user'], config['password'], config['data_server_path'], path, logger)

def _get_latest_stats_file(config):
    path = '/'.join(config['stats_file'].split('/')[:-1])
    file = config['stats_file'].split('/')[-1]
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_filepath', type=str, default="config/web.json", help='path to the config file containing the file locations')
    args = parser.parse_args()

    config = get_config(args.config_filepath, './microreact_pipeline/config/web_template.json')
    set_logging(config)
    
    _send_stats(config)
    # website_data = _aggregate_data(config)
    # _send_data(config, website_data)
