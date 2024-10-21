import pandas as pd
import argparse
import os
import sys
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
        '''input the outpath''',)
    args = parser.parse_args()
    return args.outPath

path=get_args()
if os.path.exists(path+'/report/ATAC/Table1.dapeak.csv') and os.path.exists(path+'/report/ATAC/Table2.celltypecount.csv'): 
    df1= pd.read_csv(open(path+'/report/ATAC/Table1.dapeak.csv'),encoding="utf-8",dtype=str,)
    df1 = df1.dropna()
    fw = open(path+'/report/table/peak-table.txt','w')
    for index, row in df1.iterrows():
        fw.write('<tr><td>'+row['peak.ID']+'</td>'\
                +'<td>'+row['width']+'</td>'\
                +'<td>'+row['gene.symbol']+'</td>'\
                +'<td>'+row['distance']+'</td>'\
                +'<td>'+row['p_val']+'</td>'\
                +'<td>'+row['avg_log2FC']+'</td>'\
                +'<td>'+row['pct.1']+'</td>'\
                +'<td>'+row['pct.2']+'</td>'\
                +'<td>'+row['p_val_adj']+'</td>'\
                +'<td>'+row['cluster']+'</td></tr>'+'\n'

            )

    df1 = pd.read_csv(open(path+'/report/ATAC/Table2.celltypecount.csv'),encoding="utf-8",dtype=str)
    fw = open(path+'/report/table/cell-table.txt','w')
    for index, row in df1.iterrows():
        fw.write('<tr><td>'+row['predicated.cell.type']+'</td>'\
                +'<td>'+row['number']+'</td>'\
                +'<td>'+row['ratio']+'</td>'\
            
            )
