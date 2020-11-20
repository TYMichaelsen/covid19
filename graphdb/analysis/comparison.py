import json
import argparse
import pandas as pd


def get_config(path):
    with open(path) as f:
        return json.load(f)
   
def get_data(config):
    df_metadata = pd.read_csv(config['metadata'], sep='\t')
    df_clade = pd.read_csv(config['clades'], sep='\t')
    df_infection_rate = pd.read_csv(config['infection_rate'], sep='\t')
    df_centrality = pd.read_csv(config['centrality'], sep='\t')

    df = pd.merge(df_metadata, df_clade, how='left', on='strain')
    df = df.groupby(['clade']).size().reset_index(name='cases')
    df = pd.merge(df, df_infection_rate, how='left', left_on='clade', right_on='parent')
    df = pd.merge(df, df_centrality, how='left', left_on='clade', right_on='node_id')
    
    df = df[['clade', 'infection_rate', 'cases', 'centrality']]
    out = '/srv/rbd/covid19/bisystem/staging/analysis_comparison.tsv'
    df.to_csv(out, sep='\t')

    # df = df.fillna(0)   
    # print(df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()  
    config = get_config(args.config_file)
    df = get_data(config)
    # main(graph)


