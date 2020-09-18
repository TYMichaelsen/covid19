import pandas as pd
import logging

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

def get_all_grouped_by_week():
    pass

def get_all_grouped_by_reqion():
    pass

def get_all_grouped_by_age():
    pass

    