from py2neo import Graph, Node, Relationship
import csv

host = input("Please enter the neo4j database server IP (defaults to localhost): ")
host = 'localhost' if len(host.strip(' ')) == 0 else host
password = input("Please enter the neo4j user's database password: ")
graph = Graph("bolt://{}:7687".format(host), user='neo4j',
              password=password)  # this may need to change when running from the server since there it needs the instance IP
graph.delete_all()

# Create schema
tx = graph.begin()

## Relationships
ISA = Relationship.type("ISA")

## Sex's
sex_m = Node("Sex", name="M")
sex_f = Node("Sex", name="F")
tx.create(sex_m)
tx.create(sex_f)

## Age groups
age_groups = {}
with open('../bi_system/stable_dims/age_groups.tsv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        n = Node("AgeGroup", name=row[0])
        age_groups[row[0]] = n
        tx.create(n)

## Geodata
parishes = {}
with open('../bi_system/stable_dims/parish_ses.tsv') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    next(reader, None)
    for row in reader:
        n = Node("Parish", code=row[0], name=row[1], population=int(row[2]), ghetto_area=row[3])
        parishes[row[0]] = n
        tx.create(n)

municipalities = {}
with open('../bi_system/stable_dims/municipalities.tsv') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    next(reader, None)
    for row in reader:
        n = Node("Municipality", code=row[0], name=row[1], population=int(row[4]))
        municipalities[row[0]] = n
        tx.create(n)

countries = {}
dk = Node("Country", name='Denmark', nuts_code='DK0')
countries['Denmark'] = tx.create(dk)
nuts3_regions = {}
with open('../bi_system/stable_dims/nuts3_regions.csv') as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None)
    for row in reader:
        n = Node("NUTS3_Region", code=row[0], name=row[1])
        nuts3_regions[row[0]] = n
        tx.create(n)
        tx.create(Relationship(n,"PartOf",dk))


countries = {}

## Virus Strains
strains = {}

## Medical History
risk_factors = {}
with open('../bi_system/stable_dims/risk_factors.csv') as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None)
    for row in reader:
        n = Node("RiskFactor", code=row[0], name=row[1])
        risk_factors[row[0]] = n
        tx.create(n)

print("Created graph schema")


# CREATE (A0_9:AgeGroup {name: '0-9'})
# CREATE (A10_19:AgeGroup {name: '10-19'})
# CREATE (A20_29:AgeGroup {name: '20-29'})
# CREATE (A30_39:AgeGroup {name: '30-39'})
# CREATE (A40_49:AgeGroup {name: '40-49'})
# CREATE (A40_49:AgeGroup {name: '50-59'})
# CREATE (A40_49:AgeGroup {name: '60-69'})
# CREATE (A40_49:AgeGroup {name: '70-79'})
# CREATE (A40_49:AgeGroup {name: '80-89'})
# CREATE (A40_49:AgeGroup {name: '90+'})
#

# CREATE (SM:Sex {name: 'M'})
# CREATE (SF:Sex {name: 'F'})
#
# ## Countries
# CREATE (IT:Country {name:'Italy'})
# CREATE (P1:Patient {ssi_id:'SSI-123', sex:'F', age: '42', pregnancy: 0})


def make_rel(with_node, code_field_name, lookup_dict, relation_name, rel_node_label, name_field_name=None):
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


# Create data
with open('/srv/rbd/covid19/metadata/2020-05-26-07-35_metadata.tsv') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        cv_stat = row['COVID19_Status'] if len(row['COVID19_Status']) > 0 else '0'  # error correction
        age = row['ReportAge'] if row['ReportAge']=='' else int(row['ReportAge'])
        p = Node("Person", ssi_id=row['ssi_id'], age=age, COVID19_Status=cv_stat,
                 COVID19_EndDate=row['COVID19_EndDate'], isPregnant=(row['Pregnancy'] == '1'), sequenced=(row['sequenced'] == 'Yes'))
        if (row['Pregnancy'] == '1' and row['Sex'] == 'M'):
            print('anomalous case data') # TODO extract all error checking code to a separate file
            print('SSI {}, Pregnancy {}, Sex {}'.format(row['ssi_id'],row['Pregnancy'],row['Sex']))

        tx.create(p)
        if row['Sex'] == 'F':
            tx.create(Relationship(p, "ISA", sex_f))
        elif row['Sex'] == 'M':
            tx.create(Relationship(p, "ISA", sex_m))
        # else: # make error log file
        #  print('Unrecognized Sex value: {}'.format(row['Sex']))
        ag = row['ReportAgeGrp']
        if ag in age_groups.keys():
            tx.create(Relationship(p, "InGroup", age_groups[ag]))
        # TODO report missing ag

        # Parish
        parish = make_rel(with_node=p, code_field_name='Parishcode', lookup_dict=parishes, relation_name="LivesIn",
                          rel_node_label="Parish", name_field_name="ParishName")

        # Municipality
        if parish is not None:
            muni = make_rel(with_node=parish, code_field_name='MunicipalityCode', lookup_dict=municipalities, relation_name="PartOf",
                 rel_node_label="Municipality")

        # NUTS3 Region
        if muni is not None:
            region = make_rel(with_node=muni, code_field_name='NUTS3Code', lookup_dict=nuts3_regions, relation_name="PartOf",
                 rel_node_label="NUTS3_Region", name_field_name="NUTS3Text")


        # strains
        strain_name = row['lineage']
        if len(strain_name) > 0:
            make_rel(with_node=p,code_field_name='lineage',lookup_dict=strains, relation_name="HasStrain",
                     rel_node_label="Strain")

        # Risk factors
        for field_name in risk_factors.keys():
            if row[field_name] in ['SAND','TRUE']:
                tx.create(Relationship(p, "HasRisk", risk_factors[field_name]))


        # Place of infection
        country = row['PlaceOfInfection_EN']
        if len(country) > 0:
            make_rel(with_node=p,code_field_name='PlaceOfInfection_EN',lookup_dict=countries, relation_name="PlaceOfInfection",
                     rel_node_label="Country")

# p = Node("Person", ssi_id="example_id")
# tx.create(Relationship(p,"ISA",sex_m))
# tx.create(Relationship(p,"ISA",age_groups['40-49']))
# p2 = Node("Person", ssi_id="example_id2")
# tx.create(Relationship(p2,"ISA",sex_m))
# tx.create(Relationship(p2,"ISA",age_groups['40-49']))
tx.commit()

print("Successfully created graph from metadata file.")
