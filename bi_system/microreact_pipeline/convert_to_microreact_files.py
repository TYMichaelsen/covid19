import csv
import re
import logging 

from pandas import read_excel
from datetime import date, datetime, timedelta
from enum import Enum

LOGGER = logging.getLogger("to microreact")

class FIELD():
      ID = "ID"
      orig_id = "orig_id"
      sample_date = "sample_date"
      epi_week= "epi_weel"
      country="country"
      region="region__autocolor"
      lineage="lineage__autocolor"
      latitude="latitude"
      longitude="longitude"
      day="day"
      month="month"
      year="year"

def execute_query(connection, query):
      cursor = connection.cursor()
      cursor.execute(query)
      records = cursor.fetchall()

      assert len(records) > 0, "Query resulted in 0 records."
      return records

def convert_to_microreact_format(data):
      formatted_data = []
      epidemic_start_date = _get_epi_start_date(data)
      for e in data:
            week_start_date = _get_first_day_of_week(e[2])
            data_obj = {
                  FIELD.ID:e[1],
                  FIELD.orig_id:e[0],
                  FIELD.sample_date:e[2].isocalendar()[1],
                  FIELD.epi_week:_get_epi_week(e[2], epidemic_start_date),
                  FIELD.country:"DK01",
                  FIELD.region:e[4],
                  FIELD.lineage:e[3],
                  FIELD.latitude:e[6],
                  FIELD.longitude:e[5],
                  FIELD.day:week_start_date.day,
                  FIELD.month:week_start_date.month,
                  FIELD.year:week_start_date.year
            }
            formatted_data.append(data_obj)
      return formatted_data

def save_tree(config, tree):
      path = config['out_react_nwk']
      LOGGER.info('Saving tree to {}'.format(path))
      with open(path, "w") as f:
            f.write(tree)

def get_tree(config):
      path = config['clade_tree_path']
      with open(path, "r") as f:
            return f.read()

def replace_tree_ids(data, tree):
    for e in data:
        replacement_id = e[FIELD.ID]
        original_id = e[FIELD.orig_id]

        match_idx = _find_replacement_idx(tree, original_id)
        if match_idx == -1:
            LOGGER.warning("Failed to find and replace ID: {}".format(original_id))
            continue

        tree = tree[:match_idx] + replacement_id + tree[match_idx + len(original_id):]
    return tree

def replace_data_ids(data, linelist):
      for _, row in linelist.iterrows():
            LOGGER.debug(row['gisaid_id'])

def filter_data_by_min_cases(data, config, min_cases=3):
      cases = _get_cases_per_region_week(config)
      filtered_data = []
      skipped_ids = []
      for example in data:
            match = cases.loc[(cases['Week'] == example[FIELD.sample_date]) & (cases['NUTS3Code'] == example[FIELD.region])]
            if len(match.index) != 1:
                  LOGGER.warning("Cases did not properly match the example, continuing... (cases: {}, id: {})".format(len(match.index), example[FIELD.orig_id]))
                  skipped_ids.append(example[FIELD.ID])
                  continue
            if match['Cases'].iloc[0] < min_cases:
                  skipped_ids.append(example[FIELD.ID])
                  continue
            filtered_data.append(example)
      return filtered_data, skipped_ids

def get_unmatched_ids_in_tree(tree, id_prefix_lst):
    ids = []
    keys = ["SSI-", "HH-", "Wuhan/", "AAU-"]
    for key in keys:
        for match in re.finditer(key, tree):
            start_idx = match.start()
            end_idx = tree[start_idx:start_idx + 50].index(":") + start_idx
            ids.append(tree[start_idx:end_idx])
    return ids

def add_empty_records(data, skipped_ids):
      for skipped_id in skipped_ids:
            empty_data_obj = {
                  FIELD.ID:skipped_id,
                  FIELD.orig_id:None,
                  FIELD.sample_date:None,
                  FIELD.epi_week:None,
                  FIELD.country:None,
                  FIELD.region:None,
                  FIELD.lineage:None,
                  FIELD.latitude:None,
                  FIELD.longitude:None,
                  FIELD.day:None,
                  FIELD.month:None,
                  FIELD.year:None
            }
            data.append(empty_data_obj)
      return data

def save_csv(config, data):
      path = config['out_react_tsv']
      LOGGER.info('Saving data to {}'.format(path))
      with open(path, "w") as f:
            writer = csv.DictWriter(f,
                  fieldnames=[FIELD.ID, FIELD.sample_date, FIELD.epi_week, FIELD.country, FIELD.region, FIELD.lineage, FIELD.latitude, FIELD.longitude, 
                        FIELD.day, FIELD.month, FIELD.year])
            writer.writeheader()
            for data_obj in data:
                  del data_obj[FIELD.orig_id]
                  writer.writerow(data_obj)

def get_linelist(config):
      linelist = read_excel(config['raw_ssi_file'])
      assert linelist.empty == False
      return linelist

def _get_first_day_of_week(date):
      return date - timedelta(days=date.weekday())

def _get_epi_week(infected_date, epidemic_start_date):
      w_start_infected = infected_date - timedelta(days=infected_date.weekday())
      w_start_epidemic = epidemic_start_date - timedelta(days=epidemic_start_date.weekday())
      return (_get_first_day_of_week(w_start_infected) - _get_first_day_of_week(w_start_epidemic)).days / 7

def _get_epi_start_date(data): 
      return min(e[2] for e in data)

def _datestr_to_week_func():
      return lambda date: datetime.strptime(date, '%Y-%m-%d').isocalendar()[1]
      
def _get_cases_per_region_week(config):
      linelist = get_linelist(config)
      linelist['Week']=linelist['SampleDate'].apply(_datestr_to_week_func())
      return linelist.groupby(['Week', 'NUTS3Code']).size().reset_index(name="Cases")

def _find_replacement_idx(tree, key):
      for match in re.finditer(key, tree):
           if tree[match.end()] != ":":
                continue
           return match.start()
      return -1

