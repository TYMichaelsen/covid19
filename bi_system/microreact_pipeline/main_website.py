import argparse
import logging

from pandas import read_csv
from to_website_data_files import save_website_files, upload_web_files
from config_controller import get_config
from utilities import save_df, stfp_file, datestr_to_week_and_year_func
from web_data_processing import get_all_grouped_by_age_df, get_sequenced_grouped_by_lineage_age_df, get_sequenced_grouped_by_age_df, get_all_grouped_by_region_df, get_sequenced_grouped_by_lineage_region_df, get_sequenced_grouped_by_region_df, get_sequenced_grouped_by_week_df, get_sequenced_grouped_by_lineage_week_df, get_all_grouped_by_week_df

def set_logging(config):
    logging.basicConfig(level=logging.DEBUG, filename=config['log_filepath'], filemode='w')
    formatter = logging.Formatter('%(name)-12s: %(levelname)-6s %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


def aggregate_data(config):
    logger = logging.getLogger("web data")
    df = read_csv(config['data_file'], sep="\t")

    df = get_sequenced_grouped_by_week_df(df, config)
    print(df.head())

    df = get_sequenced_grouped_by_lineage_week_df(df, config)
    print(df.head())

    df = get_all_grouped_by_week_df(df, config)
    print(df.head())

    df = get_sequenced_grouped_by_region_df(df, config)
    print(df.head())

    df = get_sequenced_grouped_by_lineage_region_df(df, config)
    print(df.head())

    df = get_all_grouped_by_region_df(df, config)
    print(df.head())

    df = get_sequenced_grouped_by_age_df(df, config)
    print(df.head())

    df = get_sequenced_grouped_by_lineage_age_df(df, config)
    print(df.head())

    df = get_all_grouped_by_age_df(df, config)
    print(df.head())

def send_stats(config):
    logger = logging.getLogger("web stats")
    df = read_csv(config['stats_file'], sep="\t")
    save_df(df, config['stats_save_path'])
    stfp_file(config['host'], config['user'], config['password'], config['stats_server_path'], config['stats_save_path'], logger)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_filepath', type=str, default="config/web.json", help='path to the config file containing the file locations')
    args = parser.parse_args()

    config = get_config(args.config_filepath, './microreact_pipeline/config/web_template.json')
    set_logging(config)
    #send_stats(config)
    aggregate_data(config)
