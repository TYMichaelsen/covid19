import math
import json
import pandas as pd
import argparse
from py2neo import Graph, Node, Relationship, NodeMatcher

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
    df_metadata.age = df_metadata.age.apply(lambda x: '' if math.isnan(x) else str(x))
    df_metadata['age_group'] = df_metadata.age.apply(lambda x: '0-9' if len(x) == 0 else f'{x[0]}0-{x[0]}9')

    #merge dfs
    df = pd.merge(df_metadata, df_clades, how='left', on='strain')
    df = pd.merge(df, df_regions, how='left', left_on=df.division, right_on=df_regions.name)
    
    #merge linelist
    ll_columns = df_linelist.columns.difference(df.columns)
    df = pd.merge(df, df_linelist[ll_columns], how='left', left_on=df.strain, right_on='ssi_id')

    #filter columns
    meta_clade_columns = ['strain', 'clade', 'virus', 'date', 'division', 'location', 'region_exposure','division_exposure', 'country_exposure', 'host', 'age', 'age_group', 'sex']
    linelist_columns = ['COVID19_Status','COVID19_EndDate','SymptomsStartDate','Pregnancy','Symptoms','Travel','ContactWithCase','Doctor','Nurse','HealthAssist','Death60Days_final','DateOfDeath_final','Occupation','CitizenshipCode','Reg_RegionCode']
    filter_columns = meta_clade_columns + linelist_columns
    df = df[filter_columns]
    
    #fix values
    df.sex = df.sex.apply(lambda x: x if type(x) == str else '')

    df.Travel = df.Travel.apply(lambda x: True if x == 'Y' else x)
    df.Travel = df.Travel.apply(lambda x: False if x == 'N' else x)

    df.ContactWithCase = df.ContactWithCase.apply(lambda x: True if x == 'Y' else x)
    df.ContactWithCase = df.ContactWithCase.apply(lambda x: False if x == 'N' else x)
    df.COVID19_Status = df.COVID19_Status.apply(lambda x: True if x ==1.0 else False)
    df.Pregnancy = df.Pregnancy.apply(lambda x: True if x == 1.0 else False)

    allowed_regions = df_regions['name'].to_list()
    df.division = df.division.apply(lambda x: x if x in allowed_regions else '')
    df.clade = df.clade.apply(lambda x: x if type(x) == str else '')
    df.clade = df.clade.apply(lambda x: x.split('/')[0])
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

def construct_nodes(g, df, df_dl):
    [g.create(Node('AgeGroup', name=v)) for v in df.age_group.unique()]
    [g.create(Node('Sex', name=v)) for v in df.sex.unique() if v != '']
    [g.create(Node('Division', name=v)) for v in df.division.unique() if v != '']
    [g.create(Node('Location', name=v)) for v in df.location.unique() if v != '']
    [g.create(Node('Strain', name=v)) for v in df.clade.unique() if v != '']
    
    g = construct_people(df, df_dl, g)
    g = construct_division_locations(df_dl, g)

    g.commit()

def construct_division_locations(df_dl, g):
    matcher = NodeMatcher(g)
    for _,row in df_dl.iterrows():
        location_node = matcher.match('Location', name=row.location).first()
        division_node = matcher.match('Division', name=row.division).first()
        if location_node == None or division_node == None:
            continue
        g.create(Relationship(location_node, 'PartOf', division_node))
    return g

def construct_people(df, df_dl, g):
    locations = df_dl.location.unique()
    matcher = NodeMatcher(g)

    for _,v in df.head(100).iterrows():
        try:
            person_node = construct_person_node(v)
            g.create(person_node)

            sex_node = matcher.match('Sex', name=v.sex).first()
            age_node = matcher.match('AgeGroup', name=v.age_group).first()
            location_node = matcher.match('Location', name=v.location).first()
            division_node = matcher.match('Division', name=v.division).first()
            strain_node = matcher.match('Strain', name=v.clade).first()

            if v.location not in locations:
                print(f'Person {v.strain} has an unregistered location {v.location}')
            else:
                g.create(Relationship(person_node, 'LivesIn', location_node))

            g.create(Relationship(person_node, 'ISA', sex_node))
            g.create(Relationship(person_node, 'InGroup', age_node))
            g.create(Relationship(person_node, 'PartOf', division_node))
            g.create(Relationship(person_node, 'HasStrain', strain_node))
        except Exception as e:
            print('Failed to construct person node.')
            print(e)
    return g

def construct_person_node(person):
    node = Node('Person')
    properties = ['strain', 'COVID19_Status', 'Pregnancy', 'SymptomsStartDate', 'date', 'Symptoms', 'Travel', 'ContactWithCase', 'Doctor', 'Nurse', 'HealthAssist', 'Death60Days_final', 'DateOfDeath_final', 'Occupation', 'CitizenshipCode', 'Reg_RegionCode']
    for prop in properties:
        if is_valid_value(person[prop]):
            node[prop] = person[prop]
    return node
    
def get_graph(config):
    con = f'bolt://{config["neo4j_server"]}:7687'
    user = 'neo4j'
    password = config['noe4j_server_password']
    
    graph = Graph(con, user=user, password=password)
    graph.delete_all()
    return graph.begin()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()
    
    config = get_config(args.config_file)
    metadata, clades, linelist, regions, division_location = load_data(config)
    df = process_data(metadata, clades, linelist, regions)
    g = get_graph(config)
    construct_nodes(g, df, division_location)
