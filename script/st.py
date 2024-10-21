import plotly as py
# from plotly.graph_objs import Scatter, Layout, Data, Scattergl
import pandas as pd
import plotly.graph_objs as go
from plotly.io import *
import os
import argparse
from scipy.interpolate import make_interp_spline
import numpy as np
import plotly.express as px

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
        '''input the outpath''',)
    parser.add_argument('--ID', type=str, help=
        '''input the ID''',)
    args = parser.parse_args()
    return args.outPath,args.ID


config = {
    'modeBarButtonsToRemove': ["autoScale2d","hoverClosestCartesian", "hoverCompareCartesian", "lasso2d",
    "zoomIn2d", "zoomOut2d", "sendDataToCloud",
    "toggleSpikelines" ,"logo"],
    'displaylogo': False,}

if __name__ == '__main__':
    path , ID= get_args()    
    df=pd.read_csv(path+"/ATAC/02.d2cfile/"+ID+".sequenceSaturation.tsv",sep="\t")

    x=df['mean_frags_per_cell']
    y=df['saturation']*100
    if len(df) > 2:
        xnew = np.linspace(x.min(),x.max(),300)
        #import statsmodels.api as sm
        #lowess = sm.nonparametric.lowess
        #ynew = lowess(y, x, frac=0.27)
        ynew = make_interp_spline(x,y)(xnew)
        #fig = px.line(df, x=ynew[:,0], y=ynew[:,1])
    else:
        xnew = x
        ynew = y
    fig = px.line(df, x=xnew, y=ynew )

    config=config
    fig.update_layout(
        autosize=False,
        width=565,
        height=500,
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray',title="Mean Fragments per Cell"),
        yaxis=dict(gridcolor='lightgray',title="Sequencing Saturation"),
        yaxis_range=[0,100]
        )
    fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
    fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
    fig.update_traces(line=dict(color="#337ab7", width=3))
    py.offline.plot(fig, filename = path+"/Joint/report/div/saturation_atac.html",auto_open=False,config=config)
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(path+"/Joint/report/div/saturation_atac.div",'w')
    fw.write(fig2)
    fw.close()
    
    # data=Scattergl(x=df["mean_frags_per_cell"],y=df["saturation"], mode='lines')
    
    # layout = Layout(xaxis=dict(
    #                     gridcolor="lightgrey",
    #                     title="mean_frags_per_cell",
    #                     color="black",
    #                     showline=True,
    #                     zeroline=True,
    #                     linewidth=1,fixedrange= True,
    #                     linecolor="black"
    #                     ),
    #                 yaxis = dict(
    #                     title="saturation",
    #                     gridcolor="lightgrey",
    #                     linewidth=1,fixedrange= True,
    #                     color="black",
    #                     linecolor="black"
    #                     ),
    #                     height=360,width=450,
    #                     plot_bgcolor='rgba(0,0,0,0)',
    #                     hovermode='closest',
    #                     paper_bgcolor='white',
    #                 title="saturation : "+str(df['saturation'].tolist()[-1])
    
                            
    #     )
                                
    # #layout=
    # fig=dict(data=data, layout=layout)
    # #fig.show()
    # config={'displayModeBar': False}
    # py.offline.plot(fig, filename =path+"/Joint/report/div/saturation_atac.html")
    # fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    # fw = open(path+"/Joint/report/div/saturation_atac.div",'w')
    # fw.write(fig2)

