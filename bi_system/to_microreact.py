import csv
from random import randint, choice

from .convert_metadata_to_mysql import get_connection
from mysql.connector import Error


# con_dict = {
#     'db_user' : 'mirop',
#     'db_password' : '6Xk7p6ErTSnb',
#     'mariadb_server' : '132.75.249.84'
# }

con_dict = {
    'db_user' : 'covid19',
    'db_password' : 'covid19Pass',
    'mariadb_server' : '172.21.232.156'
}

def select(con, query):
    cursor = con.cursor()
    cursor.execute(query)
    records = cursor.fetchall()

    print(cursor.rowcount)

    return records

con = get_connection(con_dict)

query = 'SELECT week(date), C.low_res_clade, R.name, count(*) as cases, R.longitude, R.latitude, day(date), month(date), year(date)' \
      'FROM Persons P JOIN Municipalities M ON P.MunicipalityCode=M.code JOIN NUTS3_Regions R ON M.region=R.code ' \
      'JOIN Clade_assignment C ON P.ssi_id =  C.ssi_id' \
      'GROUP BY week(date), C.low_res_clade, R.name' \
      ';' # TBD

result = select(con, query)

data = []
for idx,e in enumerate(result):
      data_obj = {
            "ID":idx+1,
            "sample_date":e[0],
            "epi_week":e[0] - 9, #when is start of the epidemic?
            "country":"DK01",
            "region":e[2],
            "lineage":e[1],
            "latitude":e[5],
            "longitude":e[4],
            "cases":e[3] if e[3] >= 3 else 3,
            # "cases":randint(1,100),
            "day":e[6],
            "month":e[7],
            "year":e[8]
      }
      data.append(data_obj)

output_filename = "./bi_system/stable_dims/microreact.tsv"

with open(output_filename, "w") as f:
      writer = csv.DictWriter(f,
            fieldnames=["ID", "sample_date", "epi_week", "country", "region", "lineage", "latitude", "longitude", "cases", "day", "month", "year"])
      writer.writeheader()
      for data_obj in data:
            writer.writerow(data_obj)
