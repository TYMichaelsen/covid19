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
        if with_node is None:
            print('Error creating relationship {} , source node is None'.format(relation_name))
            return
        if rel_node is None:
            print('Error creating relationship {} , target node is None, key was {}, dictionary contained {}'.format(relation_name, rel_node_key, lookup_dict))
            return
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


def create_dims(tx, clades_dict):
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
    tx.create(dk)
    countries['Denmark'] = dk
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
    for clade in clades_dict.keys():
        clade_details = clades_dict[clade]
        n = Node("Strain", name=clade)
        tx.create(n)
        for country_name in clade_details['countries']:
            row = {'country_name' : country_name} # just a hack to cause make_rel to work here
            make_rel(tx=tx,row=row,with_node=n,code_field_name='country_name',lookup_dict=countries,
                         relation_name="IDENTIFIED_IN", rel_node_label="Country")
        if clade_details['parent'] is not None:
            make_rel(tx=tx,row=clade_details,with_node=n,code_field_name='parent',lookup_dict=strains,
                         relation_name="EvolvedFrom", rel_node_label="Strain")


    # Medical History
    risk_factors = {}
    with open('../bi_system/stable_dims/risk_factors.tsv') as rskfile:
        reader = csv.reader(rskfile, delimiter='\t')
        next(reader, None)
        for row in reader:
            n = Node("RiskFactor", code=row[0], name=row[1])
            risk_factors[row[0]] = n
            tx.create(n)

    print("Created dimensions")
    return {'parishes': parishes, 'sex_m': sex_m, 'sex_f': sex_f, 'municipalities': municipalities,
            'nuts3_regions': nuts3_regions, 'risk_factors': risk_factors, 'strains': strains,
            'countries': countries, 'age_groups': age_groups, 'nursing_homes': {}, 'branches': {}, 'post_codes': {} }


def connect():
    host = input("Please enter the neo4j database server IP (defaults to localhost): ")
    host = 'localhost' if len(host.strip(' ')) == 0 else host
    password = input("Please enter the neo4j user's database password: ")
    graph = Graph("bolt://{}:7687".format(host), user='neo4j',
                  password=password)  # this may need to change when running from the server since there it needs the instance IP
    graph.delete_all()
    return graph


def load_data(graph, datafile, logwriter, clade_dict):
    tx = graph.begin()
    dims = create_dims(tx, clade_dict)

    # Create data
    with open(datafile) as file:
        reader = csv.DictReader(file)
        persons = {}
        i = 0
        for row in reader:
            i += 1
            age = row['ReportAge'] if row['ReportAge'] == '' else int(row['ReportAge'])
            cv_stat = row['COVID19_Status'] if len(row['COVID19_Status']) > 0 else '0'  # error correction
            p = Node("Person", ssi_id=row['ssi_id'], age=age, COVID19_Status=cv_stat,
                     COVID19_EndDate=row['COVID19_EndDate'], IsPregnant=(row['Pregnancy'] == '1')
                     , sequenced=(row['sequenced'] == '1'), SymptomsStartDate=row['SymptomsStartDate']
                     , SampleDate=row['SampleDate'], Symptoms=row['Symptoms'], Travel=(row['Travel'] == '1')
                     # ContactWithCase=(row['ContactWithCase'] == '1')
                     , Doctor=(row['Doctor'] == '1'),
                     Nurse=(row['Nurse'] == '1'), HealthAssist=(row['HealthAssist'] == '1'),
                     Death60Days_final=(row['Death60Days_final'] == '1'), DateOfDeath=row['DateOfDeath_final'],
                     Occupation=row['Occupation'], CitizenshipCode=row['CitizenshipCode'])
                     # Reg_RegionCode=row['Reg_RegionCode'])

            if row['Pregnancy'] == '1' and row['Sex'] == 'M':
                log_field_error('anomalous case data', i,
                                'SSI {}, Pregnancy {}, Sex {}'.format(row['ssi_id'], row['Pregnancy'], row['Sex']),
                                logwriter)  # TODO extract all error checking code to a separate file

            tx.create(p)
            persons[row['ssi_id']] = p
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

            # postcode
            make_rel(tx, row, with_node=p, code_field_name='ZipCodeCity', lookup_dict=dims['post_codes'],
                              relation_name="LivesIn",
                              rel_node_label="PostCode", name_field_name="zipcode_name")
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

            # strains - replaced by clade assignment file
            # make_rel(tx, row, with_node=p, code_field_name='lineage', lookup_dict=dims['strains'],
            #              relation_name="HasStrain",
            #              rel_node_label="Strain")

            # Risk factors
            for field_name in dims['risk_factors'].keys():
                if row[field_name] in ['SAND', 'TRUE']:
                    tx.create(Relationship(p, "HasRisk", dims['risk_factors'][field_name]))

            # Place of infection
            make_rel(tx, row, with_node=p, code_field_name='PlaceOfInfection_EN', lookup_dict=dims['countries'],
                         relation_name="PlaceOfInfection",
                         rel_node_label="Country")

            # Nursing home
            make_rel(tx, row, with_node=p, code_field_name='Plejehjemsnavn', lookup_dict=dims['nursing_homes'],
                         relation_name="RESIDENT_OF", rel_node_label="NursingHome")

            # Occupation branche
            make_rel(tx, row, with_node=p, code_field_name='branche1', lookup_dict=dims['branches'],
                         relation_name="OcccupationBranche", rel_node_label="Branche")

            # Occupation branche
            make_rel(tx, row, with_node=p, code_field_name='branche2', lookup_dict=dims['branches'],
                         relation_name="OcccupationBranche", rel_node_label="Branche")

            # Occupation branche
            make_rel(tx, row, with_node=p, code_field_name='branche3', lookup_dict=dims['branches'],
                         relation_name="OcccupationBranche", rel_node_label="Branche")


    # add clade info for persons in global assignment
    for clade in clade_dict.keys():
        if clade not in dims['strains'].keys():
            log_field_error('Clade',-1, 'Missing clade {} in clade dictionary'.format(clade), logwriter)
            continue
        clade_node = dims['strains'][clade]
        for ssi_id in clade_dict[clade]['cases']:
            if ssi_id in persons.keys():
                tx.create(Relationship(persons[ssi_id], "HasStrain", clade_node))
            else:
                log_field_error('Clade person',-1, 'Missing person {} in person dictionary {}'.format(ssi_id, persons.keys()), logwriter)
                

    tx.commit()
    print("loaded data")


def get_global_clades(cladefile, logwriter):
    clades = {}
    with open(cladefile, encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        i = 0
        for row in reader:
            i += 1
            if 'Strain' not in row:
                logwriter.writerow({'MessageType': 'Error', 'ErrorType': 'FATAL', 'Details': 'Could not find Strain field in row {}'.format(row)})
                continue
            ID = row['Strain'].split('/')[1].replace('ALAB-','')
            if ID.startswith('SSI') and not ID.startswith('SSI-'):
                ID = ID.replace('SSI','SSI-')
            if ID.startswith('HH') and not ID.startswith('HH-'):
                ID = ID.replace('HH','HH-')
            country = row['Strain'].split('/')[0]
            if country == 'Wuhan':
                country = 'China'
            clade = row['Clade']
            parent = None
            if '/' in clade:
                parent = '/'.join(clade.split('/')[:-1])
            # if len(row) > 2:
            #    parent = row['parent clades'] if row['parent clades'] != '--' else None
            clade_details = clades[clade] if clade in clades.keys() else {'countries': set(), 'parent': parent, 'cases': set()}
            clade_details['countries'].add(country)
            clade_details['parent'] = parent if parent is not None else clade_details['parent']
            clade_details['cases'].add(ID)
            clades[clade] = clade_details

    logwriter.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Finished parsing {} rows from {} resulting in {} clades'.format(i, cladefile, len(clades))})

    return clades

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='load a metadata file to neo4j')
    parser.add_argument('infile', type=str, help='path to the metadata file to be checked')
    parser.add_argument('globfile', type=str, help='path to the global clade file to be checked')
    parser.add_argument('errfile', type=str, help='path to the error file to be created')
    args = parser.parse_args()
    infile = check_file(args.infile)
    print('Validated infile as {}'.format(infile))
    gfile = check_file(args.globfile)
    print('Validated global clade file as {}'.format(gfile))
    errfile = check_file(args.errfile, True)
    print('Validated errfile as {}'.format(errfile))
    with open(errfile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, ['MessageType', 'Row', 'ErrorType', 'Details'])
        writer.writeheader()
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Started {}'.format(infile)})
        graph = connect()
        clades = get_global_clades(gfile, writer)
        load_data(graph, infile, writer, clades)
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Finished {}'.format(infile)})
