import pandas as pd
import logging

from utilities import datestr_to_week_func, nut3_to_nut2_func

LOGGER = logging.getLogger("to website")

def get_seq_grouped_by_week(data, field):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([field.epi_week]).size().reset_index(name="cases")
    LOGGER.debug(data_df.head())
    
def get_seq_grouped_by_region(data, field):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([field.region]).size().reset_index(name="cases")
    LOGGER.debug(data_df.head())

def get_seq_grouped_by_age(data, field):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby([field.age_group]).size().reset_index(name="cases")
    LOGGER.debug(data_df.head())

def get_all_grouped_by_week(linelist_data):
    linelist_data['Week']=linelist_data['SampleDate'].apply(datestr_to_week_func())
    return linelist_data.groupby(['Week']).size().reset_index(name="cases")

def get_all_grouped_by_region(linelist_data):
    linelist_data['NUTS2Code']=linelist_data['NUTS3Code'].apply(nut3_to_nut2_func())
    return linelist_data.groupby(['NUTS2Code']).size().reset_index(name="cases")

def get_all_grouped_by_age(linelist_data):
    return linelist_data.groupby(['SampleAgeGrp']).size().reset_index(name="cases")
    

    