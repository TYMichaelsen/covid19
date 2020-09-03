sql = 'SELECT week(date), C.low_res_clade, R.name, count(*) as cases, R.longitude, R.latitude ' \
      'FROM Persons P JOIN Municipalities M ON P.MunicipalityCode=M.code JOIN NUTS3_Regions R ON M.region=R.code ' \
      'JOIN Clade_assignment C ON P.ssi_id =  C.ssi_id' \
      'GROUP BY week(date), C.low_res_clade, R.name' \
      ';' # TBD
