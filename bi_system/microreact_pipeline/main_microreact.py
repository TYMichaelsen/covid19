import csv
import argparse
import logging
import pandas as pd

from config_controller import get_config, set_config_nextstrain, set_config_linelist, update_latest_nextstrain
from data_cleansing_metadata import check_file, check_errors
from convert_to_microreact_files import execute_query, convert_to_microreact_format, get_tree, replace_tree_ids, filter_data_by_min_cases
from convert_to_microreact_files import get_unmatched_ids_in_tree, add_empty_records, save_csv, save_tree
from upload_microreact import upload

def set_logging(config):
    logging.basicConfig(level=logging.DEBUG, filename=config['microreact_log_path'], filemode='w')
    formatter = logging.Formatter('%(name)-12s: %(levelname)-6s %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def create_metadata_files(config):
    logger = logging.getLogger("metadata")
    
    input_file = check_file(config['sequenced_metadata_file'])
    logger.info('Validated infile as {}'.format(input_file))
    
    output_file = check_file(config['staged_metadata_file'], True)
    logger.info('Validated outfile as {}'.format(output_file))

    error_file = check_file(config['cleansing_log_file'], True)
    logger.info('Validated log file as {}'.format(error_file))

    with open(error_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, ['MessageType', 'Row', 'ErrorType', 'Details'])
        writer.writeheader()
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Started {}'.format(input_file)})
        check_errors(input_file, output_file, writer)
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Finished {}'.format(input_file)})
    logger.info('Done.')

def get_data(config):
    df_metadata = pd.read_csv(config['staged_metadata_file'], sep=',')
    df_clade = pd.read_csv(config['global_clade_assignment_file'], sep='\t')
    df_regions = pd.read_csv('stable_dims/nuts2_regions.tsv', sep='\t')

    df_metadata['NUTS2_code'] = df_metadata['NUTS3Code'].apply(lambda x: str(x)[:-1])
    df = pd.merge(df_metadata, df_clade, how='inner', left_on='ssi_id', right_on='strain')
    df = pd.merge(df, df_regions, how='left', left_on='NUTS2_code', right_on='code')
    return df

def convert_to_microreact(df, tree, config):
    logger = logging.getLogger("to microreact")

    df = convert_to_microreact_format(df)
    tree = replace_tree_ids(df, tree)
    data, skipped_ids = filter_data_by_min_cases(df, config, min_cases=3)
    data = add_empty_records(data, skipped_ids)
    logger.info("Processed {}/{}".format(len(data) - len(skipped_ids), len(data)))
    return data, tree

def save_micro_react_files(config, data, tree):
    save_csv(config, data)
    save_tree(config, tree)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_filepath', type=str, default='config/microreact.json', help='path to the config file containing the file locations')
    parser.add_argument('--date', type=str, default="latest", help="nextstrain date")
    parser.add_argument('--date_folder_suffix', type=str, default="nextstrain", help="date folder suffix (e.g. 'nextstrain' in '2020-05-10_nextstrain')")
    args = parser.parse_args()

    config = get_config(args.config_filepath, './microreact_pipeline/config/microreact_template.json')
    date_str = args.date
    date_suffix = args.date_folder_suffix

    set_logging(config)
    update_latest_nextstrain(config)
    
    config = set_config_nextstrain(config, date_str, date_suffix)
    config = set_config_linelist(config, date_str)

    create_metadata_files(config)
    data = get_data(config)
    tree = get_tree(config)
    data, tree = convert_to_microreact(data, tree, config)
    
    save_micro_react_files(config, data, tree)
    upload(config)