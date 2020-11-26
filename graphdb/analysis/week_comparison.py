import json
import pandas as pd
import math
import argparse
from datetime import datetime, timedelta, date
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import QuadMesh
from matplotlib.text import Text

def get_df(path):
    df = pd.read_csv('temps/linelist.tsv', sep='\t')
    df_municipalities = pd.read_csv('../../bi_system/stable_dims/municipalities.tsv', sep='\t')
    df_municipalities = df_municipalities[['code','name']]
    df_municipalities = df_municipalities.rename(columns={'name':'MunicipalityName'})
    df = pd.merge(df, df_municipalities, how='left', left_on='MunicipalityCode', right_on='code')
    df = df.drop(columns=['code'])
    return df

def get_figure_dims(df):
    # rows = len(df)
    # cols = len(df.columns)
    # return rows / 100 * 20, cols
    return 9,15

def set_total(df):
    total_row = df.iloc[-1, :-1]
    total_col = df.iloc[:-1, -1]
    df_center = df.iloc[:-1, :-1]

    sum_x = total_row.to_numpy().sum()
    sum_y = total_col.to_numpy().sum()
    sum_c = df_center.to_numpy().sum()
    
    assert sum_x == sum_y == sum_c, 'Total cases in dim-1, dim-2 and heatmap dont match.'
    df.iloc[-1,-1] = sum_c
    return df

def get_week_interval_dates(datestr):
    date = datetime.strptime(datestr + '/1', '%G/%V/%u')
    s_date = date - timedelta(days=date.isoweekday() % 7)
    s_date = s_date + timedelta(days=1)
    e_date = s_date + timedelta(days=6)
    return s_date, e_date

def get_dims(df, attr_x, attr_y):
    values_x = df[attr_x].unique()
    values_y = df[attr_y].unique()
    values_x = [e for e in values_x if str(e).lower() != 'nan']
    values_y = [e for e in values_y if str(e).lower() != 'nan']
    
    type_x = df[attr_x].dtype
    type_y = df[attr_y].dtype
    xy_combinations = np.stack(np.meshgrid(values_x, values_y), -1).reshape(-1, 2)
    
    df = pd.DataFrame(xy_combinations)
    df = df.rename(columns={0:attr_x, 1:attr_y})
    df = df.astype({attr_x:type_x, attr_y:type_y})
    return df

def get_heatmap_interval(df):
    df = df.drop(df.tail(1).index)
    df = df.iloc[:, :-1]
    max_value = df.to_numpy().max()
    min_value = df.to_numpy().min()
    return min_value, max_value

def get_rect_df(df_1):
    df_2 = pd.DataFrame(df_1.sum(axis=1))
    df_3 = pd.DataFrame(df_1.sum(axis=0)).transpose().rename(index={0:'total'})
    df_1['total'] = df_2[0]
    df_3['total'] = math.nan
    df = df_1.append(df_3)
    df = set_total(df)
    return df

def change_backgroun_colors(df, ax):
    cols = df.shape[1]  
    rows = df.shape[0]
    quadmesh = ax.findobj(QuadMesh)[0]
    facecolors = quadmesh.get_facecolors()
    facecolors[np.arange(cols-1,cols*rows,cols)] = np.array([1,1,1,1])
    facecolors[np.arange(cols*(rows-1),cols*rows,1)] = np.array([1,1,1,1])
    quadmesh.set_facecolors = facecolors

def save_fig(fig, ax, df, title, xlbl, ylbl, output_dir):
    change_backgroun_colors(df, ax)
    change_text_color(ax)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top') 

    plt.xlabel(xlbl, fontsize=10)
    plt.ylabel(ylbl, fontsize=10)
    fig.suptitle(title, fontsize=20)
    fig.savefig(output_dir)
    plt.close()

def change_text_color(ax):
    for i in ax.findobj(Text):
        i.set_color('black')    

def main(df, xcol, xlbl, ycol, ylbl):
    df['date_linelist'] = df['date_linelist'].apply(lambda x:datetime.strptime(x, '%Y-%m-%d'))
    df['week'] = df['date_linelist'].apply(lambda x: x.strftime("%G/%V"))
    weeks = sorted(df['week'].unique())

    df_dims = get_dims(df, ycol,xcol)
    
    df_previous = None
    for week in weeks[-5:]:
        c_week_s, c_week_e = get_week_interval_dates(week)    
        df_filtered = df[(df.date_linelist >= c_week_s) & (df.date_linelist <= c_week_e)]
        
        df_grp = df_filtered.groupby([ycol, xcol]).size().reset_index(name='cases')
        df_grp = pd.merge(df_dims, df_grp, how='left', on=[ycol, xcol])
        df_grp = df_grp.fillna(0)
        df_grp = df_grp.astype({'cases':'int32'})

        df_1 = df_grp.pivot_table(values='cases',index=ycol,columns=xcol)    
        df_res = get_rect_df(df_1)

        min_value, max_value = get_heatmap_interval(df_res)
        fig_height, fig_width = get_figure_dims(df_res)

        fig, ax = plt.subplots(1, 1, figsize =(fig_width,fig_height))
        sns.heatmap(df_res, annot=True, linewidths=0.5, cmap="Greens", vmin=min_value, vmax=max_value, fmt='g')

        title = f'Danish COVID19 cases from week {c_week_s.strftime("%Y-%m-%d")} - {c_week_e.strftime("%Y-%m-%d")}'
        output_dir = f'results/fig-{xlbl}-{ylbl}-{week.split("/")[1]}.jpg'
        save_fig(fig, ax, df_res, title, xlbl, ylbl, output_dir)

        if df_previous is not None:
            df_sub = df_res.subtract(df_previous)
            min_value, max_value = get_heatmap_interval(df_sub)
            fig2, ax2 = plt.subplots(1, 1, figsize =(fig_width + 2,fig_height))
            sns.heatmap(df_sub, annot=True, linewidths=0.5, cmap="vlag", vmin=min_value, vmax=max_value, fmt='g', center=0)
            title2 = f'Danish COVID19 cases from week {c_week_s.strftime("%Y-%m-%d")} - {c_week_e.strftime("%Y-%m-%d")} compared to previous week'
            output_dir2 = f'results/fig-{xlbl}-{ylbl}-{week.split("/")[1]}-comparison.jpg'
            save_fig(fig2, ax2, df_sub, title2, xlbl, ylbl, output_dir2)
        df_previous = df_res
            



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('xcol', type=str)
    parser.add_argument('xlbl', type=str)
    parser.add_argument('ycol', type=str)
    parser.add_argument('ylbl', type=str)

    args = parser.parse_args() 

    df = get_df('temps/linelist.tsv')
    main(df, args.xcol, args.xlbl, args.ycol, args.ylbl)