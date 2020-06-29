from py2neo import Graph, Node, Relationship
import argparse
import csv
import os


def log_field_error(field_name: str, row_num: int, err_msg: str, logfilewriter: csv.DictWriter) -> None:
    """
    Utility to log parameter errors
    :param field_name: field where error was found
    :param row_num: row number where error was found
    :param err_msg
    :param logfilewriter: for output
    :return: None
    """
    logfilewriter.writerow(
        {'MessageType': 'Error', 'Row': row_num, 'ErrorType': 'Value Error',
         'Details': 'Error in {}, details: {}'.format(field_name, err_msg)})


def make_rel(tx, row, with_node, code_field_name, lookup_dict, relation_name, rel_node_label, name_field_name=None):
    rel_node = None
    rel_node_key = row[code_field_name]
    if len(rel_node_key) > 0:
        if rel_node_key in lookup_dict.keys():
            rel_node = lookup_dict[rel_node_key]
        else:
            if name_field_name is not None:
                rel_node = Node(rel_node_label, code=rel_node_key, name=row[name_field_name])
            else:
                rel_node = Node(rel_node_label, name=rel_node_key)
            lookup_dict[rel_node_key] = rel_node
            tx.create(rel_node)
        tx.create(Relationship(with_node, relation_name, rel_node))

    return rel_node


def check_file(filepath, create=False):
    """
    Checks if the supplied path contains a file, if second argument is true, creates the file if missing
    :type filepath: str
    :param filepath: relative path to check
    :type create: bool
    :param create: defaults to false , if true will create the file
    :return: absolute filepath if it does, print error message and exit otherwise
    """
    if len(filepath) == 0:
        print("invalid filepath")
        exit(-1)

    if not os.path.exists(filepath):
        if create:
            with open(filepath, 'w'):
                pass
        else:
            print("Invalid filepath: {}".format(filepath))
            exit(-1)

    return filepath


def create_dims(graph):
    # Create schema
    tx = graph.begin()

    # Relationships
    ISA = Relationship.type("ISA")

    # Sex's
    sex_m = Node("Sex", name="M")
    sex_f = Node("Sex", name="F")
    tx.create(sex_m)
    tx.create(sex_f)

    # Age groups
    age_groups = {}
    with open('../bi_system/stable_dims/age_groups.tsv') as dimfile:
        reader = csv.reader(dimfile)
        for row in reader:
            n = Node("AgeGroup", name=row[0])
            age_groups[row[0]] = n
            tx.create(n)

    # Geodata
    parishes = {}
    with open('../bi_system/stable_dims/parish_ses.tsv') as pfile:
        reader = csv.reader(pfile, delimiter='\t')
        next(reader, None)
        for row in reader:
            n = Node("Parish", code=row[0], name=row[1], population=int(row[2]), ghetto_area=row[3])
            parishes[row[0]] = n
            tx.create(n)

    municipalities = {}
    with open('../bi_system/stable_dims/municipalities.tsv') as mfile:
        reader = csv.reader(mfile, delimiter='\t')
        next(reader, None)
        for row in reader:
            n = Node("Municipality", code=row[0], name=row[1], population=int(row[4]))
            municipalities[row[0]] = n
            tx.create(n)

    countries = {}
    dk = Node("Country", name='Denmark', nuts_code='DK0')
    countries['Denmark'] = tx.create(dk)
    nuts3_regions = {}
    with open('../bi_system/stable_dims/nuts3_regions.tsv') as rfile:
        reader = csv.reader(rfile, delimiter='\t')
        next(reader, None)
        for row in reader:
            n = Node("NUTS3_Region", code=row[0], name=row[1])
            nuts3_regions[row[0]] = n
            tx.create(n)
            tx.create(Relationship(n, "PartOf", dk))

    # Virus Strains
    strains = {}

    # Medical History
    risk_factors = {}
    with open('../bi_system/stable_dims/risk_factors.tsv') as rskfile:
        reader = csv.reader(rskfile, delimiter='\t')
        next(reader, None)
        for row in reader:
            n = Node("RiskFactor", code=row[0], name=row[1])
            risk_factors[row[0]] = n
            tx.create(n)

    tx.commit()
    print("Created dimensions")
    return {'parishes': parishes, 'sex_m': sex_m, 'sex_f': sex_f, 'municipalities': municipalities,
            'nuts3_regions': nuts3_regions, 'risk_factors': risk_factors, 'strains': strains, 'countries': countries}


def load_data(datafile, logwriter):
    host = input("Please enter the neo4j database server IP (defaults to localhost): ")
    host = 'localhost' if len(host.strip(' ')) == 0 else host
    password = input("Please enter the neo4j user's database password: ")
    graph = Graph("bolt://{}:7687".format(host), user='neo4j',
                  password=password)  # this may need to change when running from the server since there it needs the instance IP
    graph.delete_all()
    dims = create_dims(graph)

    # Create data
    tx = graph.begin()
    with open(datafile) as file:
        reader = csv.DictReader(file, delimiter='\t')
        i = 0
        for row in reader:
            i += 1
            age = row['ReportAge'] if row['ReportAge'] == '' else int(row['ReportAge'])
            cv_stat = row['COVID19_Status'] if len(row['COVID19_Status']) > 0 else '0'  # error correction
            p = Node("Person", ssi_id=row['ssi_id'], age=age, COVID19_Status=cv_stat,
                     COVID19_EndDate=row['COVID19_EndDate'], IsPregnant=(row['Pregnancy'] == '1')
                     , sequenced=(row['sequenced'] == 'Yes'), SymptomsStartDate=row['SymptomsStartDate']
                     , SampleDate=row['SampleDate'])

            if row['Pregnancy'] == '1' and row['Sex'] == 'M':
                log_field_error('anomalous case data', i,
                                'SSI {}, Pregnancy {}, Sex {}'.format(row['ssi_id'], row['Pregnancy'], row['Sex']),
                                logwriter)  # TODO extract all error checking code to a separate file

            tx.create(p)
            if row['Sex'] == 'F':
                tx.create(Relationship(p, "ISA", dims['sex_f']))
            elif row['Sex'] == 'M':
                tx.create(Relationship(p, "ISA", dims['sex_m']))
            # else: # make error log file
            #  print('Unrecognized Sex value: {}'.format(row['Sex']))
            ag = row['ReportAgeGrp']
            if ag in dims['age_groups'].keys():
                tx.create(Relationship(p, "InGroup", dims['age_groups'][ag]))
            # TODO report missing ag

            # Parish
            parish = make_rel(tx, row, with_node=p, code_field_name='Parishcode', lookup_dict=dims['parishes'],
                              relation_name="LivesIn",
                              rel_node_label="Parish", name_field_name="ParishName")

            # Municipality
            if parish is not None:
                muni = make_rel(tx, row, with_node=parish, code_field_name='MunicipalityCode',
                                lookup_dict=dims['municipalities'], relation_name="PartOf",
                                rel_node_label="Municipality")

            # NUTS3 Region
            if muni is not None:
                make_rel(tx, row, with_node=muni, code_field_name='NUTS3Code', lookup_dict=dims['nuts3_regions'],
                         relation_name="PartOf",
                         rel_node_label="NUTS3_Region", name_field_name="NUTS3Text")

            # strains
            strain_name = row['lineage']
            if len(strain_name) > 0:
                make_rel(tx, row, with_node=p, code_field_name='lineage', lookup_dict=dims['strains'],
                         relation_name="HasStrain",
                         rel_node_label="Strain")

            # Risk factors
            for field_name in dims['risk_factors'].keys():
                if row[field_name] in ['SAND', 'TRUE']:
                    tx.create(Relationship(p, "HasRisk", dims['risk_factors'][field_name]))

            # Place of infection
            country = row['PlaceOfInfection_EN']
            if len(country) > 0:
                make_rel(tx, row, with_node=p, code_field_name='PlaceOfInfection_EN', lookup_dict=dims['countries'],
                         relation_name="PlaceOfInfection",
                         rel_node_label="Country")

    tx.commit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='load a metadata file to neo4j')
    parser.add_argument('infile', type=str, help='path to the file to be checked')
    parser.add_argument('errfile', type=str, help='path to the error file to be created')
    args = parser.parse_args()
    infile = check_file(args.infile)
    print('Validated infile as {}'.format(infile))
    errfile = check_file(args.errfile, True)
    print('Validated errfile as {}'.format(errfile))
    with open(errfile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, ['MessageType', 'Row', 'ErrorType', 'Details'])
        writer.writeheader()
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Started {}'.format(infile)})
        load_data(infile, writer)
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Finished {}'.format(infile)})
