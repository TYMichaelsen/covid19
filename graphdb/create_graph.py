import math
import json
import pandas as pd
import argparse
from collections import deque
from py2neo import Graph, Node, Relationship, NodeMatcher

NODES = {}
RELATIONSHIPS = deque([])
TRANSACTION = None
GRAPH = None

def get_config(path):
    with open(path) as f:
        return json.load(f)

def load_data(config):
    df_metadata = pd.read_csv(config['metadata'], sep='\t')
    df_clade = pd.read_csv(config['clades'], sep='\t')
    df_regions = pd.read_csv(config['regions'], sep='\t')
    df_linelist = pd.read_csv(config['linelist'], sep='\t')
    df_dl = pd.read_csv(config['division_location_mapping'], sep='\t')

    return df_metadata, df_clade, df_linelist, df_regions, df_dl

def process_data(df_metadata, df_clades, df_linelist, df_regions):

    #fix empty values, columns, before merging
    df_regions.name = df_regions.name.apply(lambda x: x.split('_')[1]) 
    df_metadata.division = df_metadata.division.apply(lambda x: x if type(x) == str else '')
    df_metadata.location = df_metadata.location.apply(lambda x: x if type(x) == str else '')
    df_metadata.age = df_metadata.age.apply(lambda x: '' if math.isnan(x) else str(x))
    df_metadata['age_group'] = df_metadata.age.apply(lambda x: '0-9' if len(x) == 0 else f'{x[0]}0-{x[0]}9')

    #merge dfs
    df = pd.merge(df_metadata, df_clades, how='left', on='strain')
    df = pd.merge(df, df_regions, how='left', left_on=df.division, right_on=df_regions.name)
    
    #merge linelist
    ll_columns = df_linelist.columns.difference(df.columns)
    df = pd.merge(df, df_linelist[ll_columns], how='left', left_on=df.strain, right_on='ssi_id')

    #filter columns
    meta_clade_columns = ['strain', 'clade', 'virus', 'date', 'division', 'location', 'region_exposure','division_exposure', 'country_exposure', 'host', 'age', 'age_group', 'sex', 'direct_aa_mutations']
    linelist_columns = ['COVID19_Status','COVID19_EndDate','SymptomsStartDate','Pregnancy','Symptoms','Travel','ContactWithCase','Doctor','Nurse','HealthAssist','Death60Days_final','DateOfDeath_final','Occupation','CitizenshipCode','Reg_RegionCode']
    filter_columns = meta_clade_columns + linelist_columns
    df = df[filter_columns]
    
    #fix values
    df.sex = df.sex.apply(lambda x: x if type(x) == str else '')

    df.Travel = df.Travel.apply(lambda x: True if x == 'Y' else x)
    df.Travel = df.Travel.apply(lambda x: False if x == 'N' else x)

    df.direct_aa_mutations = df.direct_aa_mutations.apply(lambda x: x if type(x) == str else '')
    df.ContactWithCase = df.ContactWithCase.apply(lambda x: True if x == 'Y' else x)
    df.ContactWithCase = df.ContactWithCase.apply(lambda x: False if x == 'N' else x)
    df.COVID19_Status = df.COVID19_Status.apply(lambda x: True if x ==1.0 else False)
    df.Pregnancy = df.Pregnancy.apply(lambda x: True if x == 1.0 else False)

    allowed_regions = df_regions['name'].to_list()
    df.division = df.division.apply(lambda x: x if x in allowed_regions else '')
    df.clade = df.clade.apply(lambda x: x if type(x) == str else '')
    return df

def is_valid_value(v):
    if type(v) == str:
        return v != ''
    if type(v) == float:
        return math.isnan(v) == False
    if type(v) == bool:
        return True
    if type(v) == int:
        return v != -1
    return False

def build_graph(config):
    global GRAPH
    
    con = f'bolt://{config["neo4j_server"]}:7687'
    user = 'neo4j'
    password = config['noe4j_server_password']
    
    GRAPH = Graph(con, user=user, password=password)
    GRAPH.delete_all()

def begin_translaction():
    global TRANSACTION
    global GRAPH
    
    TRANSACTION = GRAPH.begin()

def commit_transaction():
    global TRANSACTION
    TRANSACTION.commit() 

def add_node(label, v):
    global TRANSACTION
    global NODES

    node = Node(label, name=v)
   
    TRANSACTION.create(node)
    if label not in NODES:
        NODES[label] = {}
    NODES[label][v]=node

def add_node_obj(node):
    global TRANSACTION
    TRANSACTION.create(node)

def get_node(label, prop):
    try:
        return NODES[label][prop]
    except:
        return None

def add_relationship(node1, node2, label):
    global TRANSACTION
    global RELATIONSHIPS

    if node1 == None or node2 == None:
        print(f'\nFailed to construct relationship {label}, None node(s).')
        return

    relationship = Relationship(node1, label, node2)
    h = hash(relationship)
    if h not in RELATIONSHIPS:
        RELATIONSHIPS.appendleft(h)
        TRANSACTION.create(relationship)

def construct_nodes(df, df_dl):
    construct_metadata_nodes(df)
    construct_clades(df)
    construct_division_locations(df_dl)
    construct_people(df, df_dl)

def construct_metadata_nodes(df):
    print('Constructing age, sex, division & location nodes.')
    
    begin_translaction()
    [add_node('AgeGroup', v) for v in df.age_group.unique() if is_valid_value(v)]
    [add_node('Sex', v) for v in df.sex.unique() if is_valid_value(v)]
    [add_node('Division', v) for v in df.division.unique() if is_valid_value(v)]
    [add_node('Location', v) for v in df.location.unique() if is_valid_value(v)]
    add_node('Origin', 'Mink')
    commit_transaction()

def construct_division_locations(df_dl):
    l = len(df_dl)
    begin_translaction()
    for i,row in df_dl.iterrows():
        print(f'\rConstructing location -> division relationships {i+1}/{l}.', end='')  
        try:
            location_node = NODES['Location'][row.location]
            division_node = NODES['Division'][row.division]
            add_relationship(location_node, division_node, 'PartOf')
        except Exception as e:
            print(f'Failed to construct location -> division relationship ({row.location} -> {row.division})')
            print(e)
    print()
    commit_transaction()

def construct_clades(df):
    begin_translaction()
    ds = df.clade.apply(lambda x: x.split('/'))
    print('Constructing clade nodes.')

    [add_node('Strain', v) for v in ds.explode().unique() if is_valid_value(v)]

    l = len(ds)
    
    for i,clades in ds.iteritems():
        print(f'\rConstructing clade -> clade relationships {i+1}/{l}.', end='')
        if clades[0] == '':
            continue
        previous_node = None
        for clade in clades:
            if previous_node == None:
                previous_node = NODES['Strain'][clade]
                continue
            node = NODES['Strain'][clade]
            add_relationship(previous_node, node, 'Contains')
    print()
    commit_transaction()

def construct_people(df, df_dl):
    locations = df_dl.location.unique()

    l = len(df)
    begin_translaction()
    for i,v in df.iterrows():
        print(f'\rConstructing person node {i+1}/{l}.', end='')
        
        if i % 50 == 0:
            commit_transaction()
            begin_translaction()
        
        person_node = construct_person_node(v)
        add_node_obj(person_node)

        if v.location not in locations:
            print(f'\nPerson {v.strain} has an unregistered location {v.location}')
        else:
            add_relationship(person_node, get_node('Location', v.location), 'LivesIn')

        if 'I692V' in v.direct_aa_mutations:
            add_relationship(person_node, get_node('Origin', 'Mink'), 'InfectedBy')

        add_relationship(person_node, get_node('Sex', v.sex), 'ISA')
        add_relationship(person_node, get_node('AgeGroup', v.age_group), 'InGroup')
        add_relationship(person_node, get_node('Division', v.division), 'PartOf')
        [add_relationship(person_node, get_node('Strain', c), 'HasStrain') for c in v.clade.split('/')]

    print()
    commit_transaction()

def construct_person_node(person):
    node = Node('Person')
    properties = ['strain', 'COVID19_Status', 'Pregnancy', 'SymptomsStartDate', 'date', 'Symptoms', 'Travel', 'ContactWithCase', 'Doctor', 'Nurse', 'HealthAssist', 'Death60Days_final', 'DateOfDeath_final', 'Occupation', 'CitizenshipCode', 'Reg_RegionCode']
    for prop in properties:
        if is_valid_value(person[prop]):
            node[prop] = person[prop]
    return node
     
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()
    
    config = get_config(args.config_file)
    metadata, clades, linelist, regions, division_location = load_data(config)
    df = process_data(metadata, clades, linelist, regions)
    build_graph(config)
    construct_nodes(df, division_location)
