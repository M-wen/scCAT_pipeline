import os
import pandas as pd
import numpy as np

### to all hex
atacpath='/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/liuzhendan/pipeline/03.dsciCAT/barcode/whitelist_V1_30bp.txt'

atac = pd.read_csv(atacpath,sep="\t",header=None)
fw = open('hex_96barcode_map_true.tsv','w')
for index, row in atac.iterrows():
    a=str(row[0]).replace('A','0').replace('C','1').replace('G','2').replace('T','3')
    b=np.base_repr(int(a, 4), base=16)
    if len(b) < len(a)/2:
        b='0'+b

    fw.write(str(b)+'\t'+str(row[0])+'\n')
