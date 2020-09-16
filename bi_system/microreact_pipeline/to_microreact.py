import csv
import argparse
import logging

from config_controller import get_config, set_config_nextstrain, update_latest_nextstrain
from data_cleansing_metadata import check_file, check_errors
from convert_metadata_to_mysql import get_connection, create_schema, add_data, create_fk
from convert_to_microreact_files import QUERY, execute_query, convert_to_microreact_format, get_tree, replace_tree_ids, filter_data_by_min_cases
from convert_to_microreact_files import get_unmatched_ids_in_tree, add_empty_records, save_csv, save_tree

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

def convert_to_sql(config):
    connection = get_connection(config)
    create_schema(connection)
    add_data(connection, config['staged_metadata_file'], config['global_clade_assignment_file'])
    create_fk() 
    connection.close()

def convert_to_microreact(config):
    logger = logging.getLogger("to microreact")

    connection = get_connection(config)
    data = execute_query(connection, QUERY)
    data = convert_to_microreact_format(data)
    connection.close()

    tree = get_tree(config)
    tree = replace_tree_ids(data, tree)

    data, skipped_ids = filter_data_by_min_cases(data,config, min_cases=3)
    data = add_empty_records(data, skipped_ids)

    save_csv(config, data)
    save_tree(config, tree)
    logger.info("Processed {}/{}".format(len(data) - len(skipped_ids), len(data)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_filepath', type=str, default="config.json", help='path to the config file containing the file locations, see config.json.template in this directory')
    parser.add_argument('--date', type=str, default="latest", help="nextstrain date")
    parser.add_argument('--date_folder_suffix', type=str, default="nextstrain", help="date folder suffix (e.g. 'nextstrain' in '2020-05-10_nextstrain')")
    args = parser.parse_args()

    
    config = get_config(args.config_filepath)
    date_str = args.date
    date_suffix = args.date_folder_suffix

    set_logging(config)
    update_latest_nextstrain(config)
    config = set_config_nextstrain(config, date_str, date_suffix)
    create_metadata_files(config)
    convert_to_sql(config)
    convert_to_microreact(config)