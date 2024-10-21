import plotly as py
from plotly.graph_objs import Scatter, Layout, Data, Scattergl
import pandas as pd
import plotly.graph_objs as go
from plotly.io import *
import os
import argparse
import numpy as np
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
	'''input the outpath''',)
    args = parser.parse_args()
    return args.outPath
    
def plot_jaccard_knee_frag():
    outpath = get_args()
    df = pd.read_csv(open(outpath+"/report/ATAC/plot_input1_Bead_Barcode_Knee.csv"),encoding="utf-8",header=None,sep="\t") 
    seq_num = np.array(df[1]).tolist() 
    #import random
    #seq_num = random.sample(seq_num, 1000)
    seq_num.sort()
    seq_num = [i for i in seq_num if i > 100]
    header_num = [1,2,3,4,5,6,7,8,9,10]
    seq_num = header_num + seq_num
    ls_len = len(seq_num)
    sort_num = [i for i in range(ls_len+1)] 
    sort_num.sort(reverse=True) 
    cell_info = pd.read_csv(open(outpath+"/report/ATAC/4.cells.csv"),encoding="utf-8",header=None,sep=":") 
    #bead_number = int(cell_info[1][0].replace(",",""))
    try:
        bead_number = int(cell_info[1][0])
    except:
        bead_number = int(cell_info[1][0].replace(",",""))
    #print(bead_number)
    seq_sort = list(zip(sort_num,seq_num))
    #order change
    higher_than_point_trans = [i for i in seq_sort if int(i[1]) > bead_number]
    
    lower_than_point_trans = [i for i in seq_sort if i[1] <= bead_number]
    
    trace0_x = [i[0] for i in higher_than_point_trans]
    trace0_x.sort()  
    trace0_y = [i[1] for i in higher_than_point_trans]
    trace0_y.sort(reverse=True)
    
    trace1_x = [i[0] for i in lower_than_point_trans]
    trace1_x.sort()
    trace1_y = [i[1] for i in lower_than_point_trans]
    trace1_y.sort(reverse=True) 
    
    blue_line = list(zip(trace0_x, trace0_y))
    blue_line = [list(i) for i in blue_line]
    
    gray_line = list(zip(trace1_x, trace1_y))
    gray_line = [list(i) for i in gray_line]  
    
    trace0 = Scattergl(
        x = trace0_x,
        y = trace0_y,
        mode="lines",
        name="TRUE",
        line=dict(color="#005bac",width=5)
    )
    trace1 = Scattergl(
        x = trace1_x,
        y = trace1_y,
        mode="lines",
        name="FALSE",
        line=dict(color="gray",width=5)
        
    )
    config={'displayModeBar': False}
    #hovermode='closest',)
    data = [trace0, trace1]
    config={'displayModeBar': False}
    layout = Layout(
                        xaxis=dict(type="log", 
                        gridcolor="lightgrey",
                        title="Barcode in Rank-descending Order",
                        color="black",
                        showline=True,
                        zeroline=True,
                        linewidth=1,fixedrange= True,
                        linecolor="black"
                        ),
                        
                        yaxis = dict(
                        type="log",
                        title="Reads per Barcode",
                        gridcolor="lightgrey",
                        linewidth=1,fixedrange= True,
                        color="black",
                        linecolor="black"
                        ),
                        height=360,width=450,
                        plot_bgcolor='rgba(0,0,0,0)',
                        
                        hovermode='closest',
                        paper_bgcolor='white',
                        
                        legend=dict(
                        x=1,
                        y=1,
                        traceorder="normal",
                        font=dict(
                        family="Arial",
                        size=12,
                        color="black"
                        ),
                        bordercolor="Black",
                        borderwidth=0
                        ),
                        margin=dict(
                        l=0,
                        r=0,
                        b=0,
                        t=0,
                        pad=1
                        ),
                        font=dict(size=10)
    ) 
    
    data = [trace0, trace1]
    fig = dict(data=data, layout=layout)
    py.offline.plot(fig, filename = outpath+"/report/div/barcode_rank.html",auto_open=False,config=config)
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(outpath+"/report/div/barcode_rank.div",'w')
    fw.write(fig2)
if __name__ == '__main__':
    path = get_args()
    plot_jaccard_knee_frag()
