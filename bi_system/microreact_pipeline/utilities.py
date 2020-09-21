from datetime import datetime
from pandas import read_excel

def get_linelist(config):
      linelist = read_excel(config['raw_ssi_file'])
      assert linelist.empty == False
      return linelist

def datestr_to_week_func():
      return lambda date: datetime.strptime(date, '%Y-%m-%d').isocalendar()[1]

def nut3_to_nut2_func():
      return  lambda nut3: nut3[:-1] if type(nut3) is str else ''