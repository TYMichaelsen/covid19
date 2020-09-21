from datetime import datetime

def datestr_to_week_func():
      return lambda date: datetime.strptime(date, '%Y-%m-%d').isocalendar()[1]

def nut3_to_nut2_func():
      return  lambda nut3: nut3[:-1] if type(nut3) is str else ''