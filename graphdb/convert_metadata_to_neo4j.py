from py2neo import Graph, Node, Relationship
import csv

host = input("Please enter the neo4j database server IP (defaults to localhost): ")
host = 'localhost' if len(host.strip(' ')) == 0 else host
password = input("Please enter the neo4j user's database password: ")
graph = Graph("bolt://{}:7687".format(host),user='neo4j', password=password)  # this may need to change when running from the server since there it needs the instance IP
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
with open('../bi_system/stable_dims/age_groups.txt') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        n = Node("AgeGroup", name = row[0])
        age_groups[row[0]] = n
        tx.create(n)

## Parishes
parishes = {}

## Virus Strains
strains = {}

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

# Create data
with open('/srv/rbd/covid19/metadata/2020-05-26-07-35_metadata.tsv') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        cv_stat = row['COVID19_Status'] if len(row['COVID19_Status'])>0 else '0' # error correction
        p = Node("Person", ssi_id=row['ssi_id'], age=row['ReportAge'], COVID19_Status=cv_stat,
                 COVID19_EndDate=row['COVID19_EndDate'], IsPregnant=(row['Pregnancy']=='1')  )
        tx.create(p)
        if row['Sex'] == 'F':
            tx.create(Relationship(p,"ISA",sex_f))
        elif row['Sex'] == 'M':
            tx.create(Relationship(p,"ISA",sex_m))
        # else: # make error log file
            #  print('Unrecognized Sex value: {}'.format(row['Sex']))
        ag = row['ReportAgeGrp']
        if ag in age_groups.keys():
            tx.create(Relationship(p,"InGroup",age_groups[ag]))
        # TODO report missing ag

        # Parish
        parish_code = row['Parishcode']
        if len(parish_code) > 0:
            if parish_code in parishes.keys():
                parish = parishes[parish_code]
            else:
                parish = Node("Parish", code=parish_code, name=row['ParishName'])
                parishes[parish_code] = parish
                tx.create(parish)

            tx.create(Relationship(p,"LivesIn",parish))

        # strains
        strain_name = row['lineage']
        if len(strain_name) > 0:
            if strain_name in strains.keys():
                strain = strains[strain_name]
            else:
                strain = Node("Strain", name=strain_name)
                strains[strain_name] = strain
                tx.create(strain)

            tx.create(Relationship(p,"HasStrain",strain))


# p = Node("Person", ssi_id="example_id")
# tx.create(Relationship(p,"ISA",sex_m))
# tx.create(Relationship(p,"ISA",age_groups['40-49']))
# p2 = Node("Person", ssi_id="example_id2")
# tx.create(Relationship(p2,"ISA",sex_m))
# tx.create(Relationship(p2,"ISA",age_groups['40-49']))
tx.commit()

print("Successfully created graph from metadata file.")
