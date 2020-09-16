import json
from datetime import date
import csv
import os
import logging

import mysql.connector
from mysql.connector import errorcode
import argparse

LOGGER = logging.getLogger("to mysql")

DB_NAME = 'covid19'
DIM_PATH = 'stable_dims'
TABLES = {'Persons': ("CREATE TABLE `Persons` ("
                      "  `ssi_id` varchar(50) NOT NULL PRIMARY KEY,"
                      "  `age` int(11),"
                      "  `age_group` varchar(16),"
                      "  `sex` enum('M','F'),"
                      "  `COVID19_Status` enum('0','1','2') NOT NULL,"
                      "  `Parishcode` INTEGER,"
                      "  `MunicipalityCode` INTEGER,"
                      "  `lineage` VARCHAR(10),"
                      , '', ['ssi_id', 'age', 'age_group', 'sex', 'COVID19_Status', 'Parishcode', 'lineage']),
          'Countries': (
              "CREATE TABLE `Countries` (`country` varchar(35) PRIMARY KEY, population INTEGER COMMENT '(2020)',"
              " `land_area` INTEGER COMMENT '(Km²)', `density` INTEGER COMMENT '(P/Km²)')",
              'countries.tsv', ['country', 'population', 'land_area', 'density']),
          'Municipalities': ("CREATE TABLE Municipalities (code INTEGER PRIMARY KEY, name VARCHAR (35), "
                             "administrative_center VARCHAR(40), area FLOAT, population INTEGER COMMENT '(2012-01-01)', "
                             "region CHAR(5))", "municipalities.tsv", ['code', 'name', 'administrative_center',
                                                                       'area', 'population', 'region']),
          'AgeGroups': ("CREATE TABLE AgeGroups (age_group VARCHAR(5) PRIMARY KEY, meta_group VARCHAR(8))",
                        "age_groups.tsv", ['age_group', 'meta_group']),
          'NUTS3_Regions': (
              "CREATE TABLE NUTS3_Regions(code CHAR(5) PRIMARY KEY, `name` VARCHAR(20), `latitude` varchar(15), `longitude` varchar(15) )", "nuts3_regions.tsv",
              ['code', 'name', 'latitude', 'longitude' ]),
          'Parishes': (
              "CREATE TABLE Parishes(code INTEGER PRIMARY KEY, `name` VARCHAR(35))", "parish.tsv", ['code', 'name']),
          'Clade_assignment': ("CREATE TABLE Clade_assignment(`ssi_id` varchar(50) NOT NULL PRIMARY KEY,"
                               "`low_res_clade` varchar(45),"
                               "`high_res_clade`  varchar(100))",'', ['ssi_id','low_res_clade','high_res_clade'])}

BOOL_FIELDS = ['Diabet', 'Neuro', 'Cancer', 'Adipos', 'Nyre', 'Haem_c', 'Card_dis', 'Resp_dis', 'Immu_dis',
               'Other_risk', 'Pregnancy', 'Doctor', 'Nurse']
DATE_FIELDS = ['date', 'COVID19_EndDate', 'SymptomsStartDate']

def get_connection(config_dict):
    try:
        cnxn = mysql.connector.connect(user=config_dict['db_user'], password=config_dict['db_password'],
                                       host=config_dict['mariadb_server'],
                                       database=DB_NAME)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            LOGGER.error("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            LOGGER.error("Database does not exist")
        else:
            LOGGER.error(err)
    else:
        return cnxn

def check_file(filepath):
    """
    Checks if the supplied path contains a file
    :type filepath: str
    :param filepath: relative path to check
    :return: absolute filepath if it does, print error message and exit otherwise
    """
    if len(filepath) == 0 or not os.path.exists(filepath):
        LOGGER.error("file not found")
        exit(-1)
    else:
        return filepath

def create_schema(cnxn):
    cursor = cnxn.cursor()
    try:
        cursor.execute("USE {}".format(DB_NAME))
        for table_name in TABLES:
            table_description = TABLES[table_name][0]
            try:
                # TODO: drop them all like this: SELECT concat('DROP TABLE IF EXISTS `', table_name, '`;')
                # FROM information_schema.tables
                # WHERE table_schema = 'MyDatabaseName';
                LOGGER.info("Creating table {}: ".format(table_name))
                cursor.execute("DROP TABLE IF EXISTS {}".format(table_name))
                if table_name == 'Persons':
                    for df in DATE_FIELDS:
                        table_description += "{} date,".format(df)
                    for df in BOOL_FIELDS:
                        table_description += "{} varchar(3),".format(df)
                    table_description = table_description[:-1]
                    table_description += ") ENGINE=InnoDB"
                cursor.execute(table_description)
            except mysql.connector.Error as err:
                if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
                    LOGGER.error("already exists.")
                else:
                    LOGGER.error(err.msg)
            else:
                LOGGER.info("Created the schema")
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_BAD_DB_ERROR:
            LOGGER.error("Database {} does not exists.".format(DB_NAME))
            exit(-1)

def add_data(cnxn, filepath, clade_filepath):
    cursor = cnxn.cursor()

    # Dimensions
    for table_name, definition in TABLES.items():
        if table_name == 'Persons' or table_name == 'Clade_assignment':
            continue
        _, dim_filepath, field_list = definition
        field_list_str = ','.join(field_list)
        # noinspection PyUnusedLocal
        place_holders = ','.join(['%s' for fl in field_list])
        insert_string = ("INSERT INTO {} ({}) VALUES ({})".format(table_name, field_list_str, place_holders))
        with open(os.path.join(DIM_PATH, dim_filepath)) as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            for row in reader:
                data = [row[col] for col in field_list]
                try:
                    cursor.execute(insert_string, data)
                except mysql.connector.Error as err:
                    LOGGER.error(err)
                    LOGGER.error("Failed data: {} for insert statement: {}".format(data, insert_string))
        LOGGER.info("Finished loading {}".format(table_name))

    # Persons
    bool_field_names = ','.join(BOOL_FIELDS)
    bf_ss = ', '.join(['%s' for fl in BOOL_FIELDS])
    date_field_names = ','.join(DATE_FIELDS)
    df_ss = ', '.join(['%s' for fl in DATE_FIELDS])
    add_person = ("INSERT INTO Persons "
                  "(ssi_id, age, age_group, sex, COVID19_Status, Parishcode, MunicipalityCode, lineage, {}, {}) "
                  "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, {}, {})".format(bool_field_names, date_field_names, bf_ss,
                                                                       df_ss))

    people = set()
    with open(filepath) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            cv_stat = row['COVID19_Status'] if len(row['COVID19_Status']) > 0 else '0'  # error correction
            #    print("stat: {} {}".format(row['COVID19_Status'], cv_stat))
            age = None if row['ReportAge'] == '' else int(row['ReportAge'])
            ag = row['ReportAgeGrp']
            sex = row['Sex'] if row['Sex'] in ['M', 'F'] else None
            pc = int(row['Parishcode']) if len(row['Parishcode']) > 0 else None
            mc = int(row['MunicipalityCode']) if len(row['MunicipalityCode']) > 0 else None
            data_person = [row['ssi_id'], age, ag, sex, cv_stat, pc, mc, row['lineage']]
            boolean_field_data = [row[f] if len(row[f]) > 0 else None for f in BOOL_FIELDS]
            date_field_data = [get_date(row, df) for df in DATE_FIELDS]
            data_person.extend(boolean_field_data)
            data_person.extend(date_field_data)
            data_person = tuple(data_person)
            people.add(data_person[0])

            # COVID19_EndDate=row['COVID19_EndDate'], isPregnant=(row['Pregnancy'] == '1'), sequenced=(row['sequenced'] == 'Yes'))
            try:
                cursor.execute(add_person, data_person)
            except mysql.connector.Error as err:
                LOGGER.error(err)
                LOGGER.error("Failed data: {}".format(data_person))
                LOGGER.error("Insert statement:")
                LOGGER.error(add_person)

    # Clade assignment
    add_clade =  ("INSERT INTO Clade_assignment(ssi_id, low_res_clade, high_res_clade) VALUES (%s, %s, %s)")
    with open(clade_filepath) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            pid = row['strain']
            if pid is None or pid not in people:
                continue

            high_res =  row['clade'].strip(' ')
            cl_arr = high_res.split('/')
            if len(cl_arr)<=1:
                low_res = high_res
            else:
                low_res = cl_arr[0]

            data_clade = (pid, low_res, row['clade'].strip(' '))
            try:
                cursor.execute(add_clade, data_clade)
            except mysql.connector.Error as err:
                LOGGER.error(err)
                LOGGER.error("Failed data: {}".format(data_clade))
                LOGGER.error("Insert statement:")
                LOGGER.error(add_clade)



    cnxn.commit()
    LOGGER.info("Loaded all data from {}".format(filepath))

def get_date(row, field_name):
    enddate_arr = row[field_name].split('-')  # e.g. '2020-03-22'.split('-')
    try:
        if len(enddate_arr) == 3 and len(enddate_arr[0]) == 4:
            enddate = date(year=int(enddate_arr[0]), month=int(enddate_arr[1]), day=int(enddate_arr[2]))
        elif len(enddate_arr) == 3 and len(enddate_arr[2]) == 4:
            enddate = date(year=int(enddate_arr[2]), month=int(enddate_arr[1]), day=int(enddate_arr[0]))
        else:
            enddate = None
    except ValueError as err:
        LOGGER.error(err)
        LOGGER.error("Failed data: {} year {} month {} day {} ".format(enddate_arr, int(enddate_arr[0]),
                                                                int(enddate_arr[1]), int(enddate_arr[2])))
        enddate = None
    return enddate

def create_fk():
    # TODO
    LOGGER.warning("Warning: Foreign keys were not created (Not implemented).")