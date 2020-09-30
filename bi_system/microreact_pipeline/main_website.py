import argparse
import logging

from to_website_data_files import save_website_files, upload_web_files
from website_stats import send_stats
from config_controller import get_config

def set_logging(config):
    logging.basicConfig(level=logging.DEBUG, filename=config['microreact_log_path'], filemode='w')
    formatter = logging.Formatter('%(name)-12s: %(levelname)-6s %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_filepath', type=str, default="config/web.json", help='path to the config file containing the file locations')
    args = parser.parse_args()

    config = get_config(args.config_filepath, './microreact_pipeline/config/web_template.json')
    send_stats(config)
