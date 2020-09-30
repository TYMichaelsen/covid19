import logging 
import csv

LOGGER = logging.getLogger('website stat files')

def send(config):
  _load_stats_file(config['stats_file'])

def _load_stats_file(path):
  with open(path, 'r') as f:
    data = csv.reader(f, delimiter="\t")
  print(data)