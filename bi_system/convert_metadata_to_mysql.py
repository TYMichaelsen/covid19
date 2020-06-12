from datetime import date
import csv

import mysql.connector
from mysql.connector import errorcode
import argparse

DB_NAME = 'covid19'
TABLES = {'Persons': "CREATE TABLE `Persons` ("
                      "  `ssi_id` varchar(14) NOT NULL,"
                      "  `age` int(11),"
                      "  `age_group` varchar(16),"
                      "  `sex` enum('M','F'),"
                      "  `COVID19_Status` enum('0','1','2') NOT NULL,"
                      "  `COVID19_EndDate` date,"
                      "  PRIMARY KEY (`ssi_id`)"
                      ") ENGINE=InnoDB"}


def get_connection():
    try:
        host = input("Please enter the mysql database server IP (defaults to localhost): ")
        host = 'localhost' if len(host.strip(' ')) == 0 else host
        password = input("Please enter the mysql user's database password: ")
        cnxn = mysql.connector.connect(user='covid19', password=password,
                                       host=host,
                                       database=DB_NAME)

    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    else:
        return cnxn


def check_file(filepath):
    """
    Checks if the supplied path contains a file
    :type filepath: str
    :param filepath: relative path to check
    :return: absolute filepath if it does, print error message and exit otherwise
    """
    if len(filepath) == 0:
        print("file not found")
        exit(-1)
    else:
        return filepath


def create_schema(cnxn):
    cursor = cnxn.cursor()
    try:
        cursor.execute("USE {}".format(DB_NAME))
        for table_name in TABLES:
            table_description = TABLES[table_name]
            try:
                # TODO: drop them all like this: SELECT concat('DROP TABLE IF EXISTS `', table_name, '`;')
                # FROM information_schema.tables
                # WHERE table_schema = 'MyDatabaseName';
                print("Creating table {}: ".format(table_name), end='')
                cursor.execute("DROP TABLE IF EXISTS {}".format(table_name))
                cursor.execute(table_description)
            except mysql.connector.Error as err:
                if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
                    print("already exists.")
                else:
                    print(err.msg)
            else:
                print("Created the schema")
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database {} does not exists.".format(DB_NAME))
            exit(-1)




def clear_data(cnxn):
    cursor = cnxn.cursor()
    cursor.execute("USE {}".format(DB_NAME))
    for table in TABLES.keys():
        cursor.execute("TRUNCATE TABLE {};".format(table))
    cnxn.commit()
    print("Cleared existing data on MariaDB server")


def add_data(cnxn, filepath):
    cursor = cnxn.cursor()
    add_person = ("INSERT INTO Persons "
                   "(ssi_id, age, age_group, sex, COVID19_Status, COVID19_EndDate) "
                   "VALUES (%s, %s, %s, %s, %s, %s)")


    with open('/srv/rbd/covid19/metadata/2020-05-26-07-35_metadata.tsv') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            cv_stat = int(row['COVID19_Status']) if len(row['COVID19_Status']) > 0 else 0  # error correction
            age = None if row['ReportAge'] == '' else int(row['ReportAge'])
            ag = row['ReportAgeGrp']
            enddate = row['COVID19_EndDate'].split('-') # e.g. '2020-03-22'.split('-')
            enddate = date(int(enddate[0]),int(enddate[1]),int(enddate[2])) if len(enddate) == 3 else None
            data_person = (row['ssi_id'], age, ag , row['Sex'], cv_stat, enddate)
            # COVID19_EndDate=row['COVID19_EndDate'], isPregnant=(row['Pregnancy'] == '1'), sequenced=(row['sequenced'] == 'Yes'))
            cursor.execute(add_person, data_person)

    cnxn.commit()
    print("Loaded all data from {}".format(filepath))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='convert a metadata file to SQL relational format and inject to the MariaDB')
    parser.add_argument('file', type=str,
                        help='path to the file to be loaded')
    parser.add_argument('--recreate_schema', dest='re_schema', type=str, default='no',
                        help='yes to recreate the schema')

    args = parser.parse_args()
    file = check_file(args.file)
    cnx = get_connection()
    if 'y' in args.re_schema.lower():
        create_schema(cnx)

    clear_data(cnx)
    add_data(cnx, file)

    cnx.close()
