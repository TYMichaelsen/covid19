import pandas as pd
import logging

LOGGER = logging.getLogger("to website")

def get_all_grouped_by_week():
    pass

def get_seq_grouped_byt_week(data):
    data_df = pd.DataFrame(data)
    data_df = data_df.groupby(['epi_weel']).size().reset_index(name="Cases")
    LOGGER.debug(data_df.head())
    
def get_all_grouped_by_reqion():
    pass


    