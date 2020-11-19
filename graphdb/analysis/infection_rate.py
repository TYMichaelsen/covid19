import json
import argparse
import pandas as pd

from py2neo import Graph, Node, Relationship, NodeMatcher

def get_config(path):
    with open(path) as f:
        return json.load(f)

def get_graph(config):
    con = f'bolt://{config["neo4j_server"]}:7687'
    user = 'neo4j'
    password = config['noe4j_server_password']  
    return Graph(con, user=user, password=password)

def get_strain_connections(graph):
    result = graph.run("MATCH (s1:Strain)-[]->(s2:Strain) RETURN s1,s2").data()
    parents = []
    for v in result:
        new_child = v['s1']
        new_parent = v['s2']

        existing_parent = None
        for parent in parents:
            if parent.get_node() == new_parent:
                existing_parent = parent
                break
        
        if existing_parent == None:
            parent = Parent(new_parent)
            parent.add_child(new_child)
            parents.append(parent)
        else:
            existing_parent.add_child(new_child)
    return parents

def main(graph):
    # GET ALL NODES SUCH AS: PARENT ---> CHILD1, CHILD2, ...
    parents = get_strain_connections(graph)
    row_list = []
    for parent in parents:
        # GET PERSON CONNECTIONS TO PARENT (X)
        key = 'v'
        result = graph.run(f'MATCH (:Person)-[r:HasStrain]->(:Strain {{name:\"{parent.name()}\"}}) return count(r) as {key}').data()
        x = result[0][key]

        #GET PERSON CONNECTION TO EACH CHILD Y = CH1 + CH2 +...
        y = 0
        strains = []
        for child in parent.get_children():        
            strain = child['name']
            result = graph.run(f'MATCH (:Person)-[r:HasStrain]->(:Strain {{name:\"{strain}\"}}) return count(r) as {key}').data()
            y += result[0][key]
            strains.append(strain)
        
        # IR(X) = X + Y / X
        ir = x + y / x
        row = {
            'parent':parent.name(),
            'infection_rate':ir,
            'children':strains
        }
        row_list.append(row)
        # print(f'IR({parent.name()}) = {ir} {10*" "} children: {strains}')
    df = pd.DataFrame(row_list)
    df.to_csv('/srv/rbd/covid19/bisystem/staging/infection_rate.tsv', sep='\t')
    print(df)
    
class Parent():
    def __init__(self, node):
        self._node = node
        self._children = []

    def add_child(self, child_node):
        for child in self._children:
            if child['name'] == child_node['name']:
                return
        self._children.append(child_node)
    
    def name(self):
        return self.get_node()['name']

    def get_node(self):
        return self._node

    def get_children(self):
        return self._children
    
    def print(self):
        print(f"Parent: {self._node['name']}")
        for child in self._children:
            print(f"     {self._node['name']} ----> {child['name']}")
        # print()    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()  
    config = get_config(args.config_file)
    graph = get_graph(config)
    main(graph)


