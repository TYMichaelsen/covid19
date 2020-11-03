import re
import logging 
import uuid
import pandas as pd
import math

from datetime import datetime, timedelta
from utilities import datestr_to_week_func, nut3_to_nut2_func, get_linelist

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
      age_group="age_group"

def execute_query(connection, query):
      cursor = connection.cursor()
      cursor.execute(query)
      records = cursor.fetchall()

      assert len(records) > 0, "Query resulted in 0 records."
      return records

def convert_to_microreact_format(df):
      LOGGER.info('Converting metadatada to microreact format...')
      df['date'] = df['date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
      df[FIELD.ID] = pd.Series(index=range(len(df.index))).apply(lambda _: str(uuid.uuid4()))
      df[FIELD.sample_date] = df['date'].apply(lambda x: x.isocalendar()[1])
      df[FIELD.epi_week] = df['date'].apply(_get_epi_week, args=(min(df['date']),))
      df[FIELD.country] = pd.Series(index=range(len(df.index))).apply(lambda _: 'DK01')
      df[FIELD.day] = df['date'].apply(lambda x: _get_first_day_of_week(x).day)
      df[FIELD.month] = df['date'].apply(lambda x: _get_first_day_of_week(x).month)
      df[FIELD.year] = df['date'].apply(lambda x: _get_first_day_of_week(x).year)
      df[FIELD.age_group] = df['ReportAgeGrp'].astype('str')
      df[FIELD.lineage] = df['clade'].apply(lambda x: x.split('/')[0])

      df = df.rename(columns={'ssi_id':FIELD.orig_id, 'NUTS2_code':FIELD.region})
      df = df[[FIELD.ID, FIELD.orig_id, FIELD.sample_date, FIELD.epi_week, FIELD.country, FIELD.region, FIELD.lineage, FIELD.latitude, FIELD.longitude, FIELD.day, FIELD.month, FIELD.year, FIELD.age_group]]  
      return df

def save_tree(config, tree):
      path = config['out_react_nwk']
      LOGGER.info('Saving tree to {}'.format(path))
      with open(path, "w") as f:
            f.write(tree)

def get_tree(config):
      path = config['clade_tree_path']
      with open(path, "r") as f:
            return f.read()

def replace_tree_ids(df, tree):
      LOGGER.info('Replacing Tree IDs with generrated UUIDs...')
      for _,row in df.iterrows():
            replacement_id = row[FIELD.ID]
            original_id = row[FIELD.orig_id]

            match_idx = _find_replacement_idx(tree, original_id)
            if match_idx == -1:
                  LOGGER.warning("Failed to find and replace ID: {}".format(original_id))
                  continue

            tree = tree[:match_idx] + replacement_id + tree[match_idx + len(original_id):]
      return tree

def filter_data_by_min_cases(df, config, min_cases=3):
      LOGGER.info('Filtering data with minimum number of cases per week and region (min number of cases - {})...'.format(min_cases))
      cases = _get_cases_per_region_week(config)
      filtered_data = []
      skipped_ids = {}

      region_less = []

      for _,example in df.iterrows():
            if example[FIELD.region] == 'Unknown':
                  region_less.append(example.to_dict())
                  continue
            
            match = cases.loc[(cases[FIELD.sample_date] == example[FIELD.sample_date]) & (cases[FIELD.region] == example[FIELD.region])]
            if len(match.index) != 1:
                  LOGGER.warning("Cases did not properly match the example, continuing... (cases: {}, id: {})".format(len(match.index), example[FIELD.orig_id]))
                  skipped_ids[example[FIELD.ID]] = example[FIELD.lineage]
                  continue
            
            if match['Cases'].iloc[0] < min_cases:
                  LOGGER.warning('Example with id: {} does not satisfy minimum number of cases. Cases: {}'.format(example[FIELD.orig_id], len(match.index)))
                  skipped_ids[example[FIELD.ID]] = example[FIELD.lineage]
                  continue

            filtered_data.append(example.to_dict())

      if len(region_less) < min_cases:
            LOGGER.warning('There are less than {} cases without region - skipping...'.format(min_cases))
            for example in region_less:
                  skipped_ids[example[FIELD.ID]] = example[FIELD.lineage]
      else:
            LOGGER.info('Appending region less cases...')
            filtered_data.extend(region_less)

      LOGGER.info('Data entries after filtering: {}'.format(len(filtered_data)))
      LOGGER.info('Skipped entries: {}'.format(len(skipped_ids)))
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
                  FIELD.sample_date:'Unknown',
                  FIELD.epi_week:'Unknown',
                  FIELD.country:'DK01',
                  FIELD.region:'Unknown',
                  FIELD.lineage:skipped_ids[skipped_id],
                  FIELD.latitude:None,
                  FIELD.longitude:None,
                  FIELD.day:None,
                  FIELD.month:None,
                  FIELD.year:None,
                  FIELD.age_group:None
            }
            data.append(empty_data_obj)
      return data

def save_csv(config, data):
      path = config['out_react_tsv']
      LOGGER.info('Saving data to {}'.format(path))

      df = pd.DataFrame(data)
      df = df.drop(columns=[FIELD.orig_id, FIELD.age_group], axis=1)
      df.to_csv(path, sep='\t')

def _get_first_day_of_week(date):
      return date - timedelta(days=date.weekday())

def _get_epi_week(infected_date, epidemic_start_date):
      w_start_infected = infected_date - timedelta(days=infected_date.weekday())
      w_start_epidemic = epidemic_start_date - timedelta(days=epidemic_start_date.weekday())
      return (_get_first_day_of_week(w_start_infected) - _get_first_day_of_week(w_start_epidemic)).days / 7

# def _get_epi_start_date(df): 
#       return min(df['date'])
#       # return min(e[2] for e in data)

def _get_cases_per_region_week(config):
      df = get_linelist(config)
      df[FIELD.sample_date] = df['SampleDate'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
      df[FIELD.sample_date] = df[FIELD.sample_date].apply(lambda x: x.isocalendar()[1])
      df[FIELD.region] = df['NUTS3Code'].apply(lambda x: str(x)[:-1] if 'DK' in str(x) else 'Unknown')
      return df.groupby([FIELD.sample_date, FIELD.region]).size().reset_index(name="Cases")

def _find_replacement_idx(tree, key):
      for match in re.finditer(key, tree):
           if tree[match.end()] != ":":
                continue
           return match.start()
      return -1

