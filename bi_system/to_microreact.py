sql = 'SELECT week(date), lineage, R.name, count(*) as cases, R.longitude, R.latitude ' \
      'FROM Persons P JOIN Municipalities M ON P.MunicipalityCode=M.code JOIN NUTS3_Regions R ON M.region=R.code ' \
      'GROUP BY week(date), lineage, R.name' \
      ';' # TBD
