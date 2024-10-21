import pandas as pd
import os
import argparse
import re
import datetime

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
	'''input the outpath''',)
    parser.add_argument('--htmlTemplate', type=str, help=
	'''input the html template''',)
    parser.add_argument('--ID', type=str, help=
        '''input the ID''',)
    parser.add_argument('--intron', type=str, help=
        '''True or False of intron reads''',)
    args = parser.parse_args()
    return [args.outPath, args.htmlTemplate, args.ID, args.intron]
 
def get_args_from_file():
    path=get_args()[0]
    intron=get_args()[3]
    csv = [path+'/report/RNA/1.cell_report.csv',\
    path+'/report/RNA/3.cDNA_sequencing_report.csv',\
    path+'/report/RNA/4.alignment_report.csv',\
    path+'/report/RNA/5.anno_report.csv',\
    path+'/report/ATAC/1.cell_report.csv',\
    path+'/report/ATAC/2.sample.csv',\
    path+'/report/ATAC/3.sequencing.csv',\
    path+'/report/ATAC/4.cells.csv',\
    path+'/report/ATAC/5.library.QC.csv',path+'/report/ATAC/5_2.library.QC.csv',path+'/report/ATAC/3.mapping.csv',\
    path+'/report/RNA/2.sample.csv',path+'/report/saturation.csv'
    ]
    
    stat = dict()
    if intron=='true':
        stat['intron_boolean'] = 'True'
    else:
        stat['intron_boolean'] = 'False'
    for i in range(len(csv)):
        if i==0:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['samplename'] = df[1][0]
            stat['species'] = df[1][1]
            stat['estm_Num_cell'] = df[1][2]
            stat['mean_r_per_c'] = format(int(df[1][3]),",")
            stat['mean_UMI_per_c'] = format(int(df[1][4]),",")
            stat['median_UMI_per_c'] = format(int(df[1][5]),",")
            stat['total_gene'] = format(int(df[1][6]),",")
            stat['mean_genes_per_c'] = format(int(df[1][7]),",")
            stat['median_genes_per_c'] = format(int(df[1][8]),",")
            stat['cluster_cell'] = format(int(df[1][9]),",")
        if i==1:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['cDNA_num_frag'] = format(int(df[1][0]),",")
            stat['cDNA_frag_pass_QC'] = df[1][1]
            stat['cDNA_frag_low_qual'] = df[1][2]
            stat['cDNA_frag_fail_bar'] = df[1][3]
            stat['cDNA_adapter_filter'] = df[1][4]
            stat['cDNA_frag_exact_bar'] = df[1][5]
            stat['cDNA_Q30_c_bar'] = df[1][7]
            stat['cDNA_Q30_s_bar'] = df[1][8]
            stat['cDNA_Q30_UMI'] = df[1][9]
            stat['cDNA_Q30_r'] = df[1][10]
        if i==2:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")        
            stat['raw_r'] = format(int(df[1][0]),",")
            stat['map_r'] = df[1][1]
            stat['plus_strd'] = df[1][2]
            stat['minus_strd'] = df[1][3]
            stat['mito_ratio'] = df[1][4]
            stat['map_qual_corrt_r'] = df[1][5]            
        if i==3:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=",")
            stat['r_m_geno'] = df[1][0]
            stat['r_m_ex'] = df[1][1]
            stat['r_m_intro'] = df[1][2]
            stat['r_m_anti'] = df[1][4]
            stat['r_m_inter'] = df[1][5]
        if i==4:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['estimate_num_of_cell'] = df[1][0]
            stat['median_frag_per_cell'] = df[1][1]
            stat['median_frac_peaks'] = str(round(float(df[1][2].strip('%')),2)) +'%'
            stat['median_frac_tss'] = str(round(float(df[1][3].strip('%')),2)) +'%'
        if i==5:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['sample_id'] = df[1][0]
            stat['fastq_path'] = df[1][1]
            stat['pipeline_version'] = df[1][2]
            stat['ref_path'] = df[1][3]
        if i==6:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")        
            stat['read_pairs'] = format(int(df[1][0]),",")
            stat['frac_valid_barcode'] = df[1][1]
        if i==7:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['bead_thres'] = format(int(df[1][0]),",")
            stat['bead_number'] = format(int(df[1][1]),",")
            stat['jaccard_thres'] = df[1][2]
        if i==8:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['frac_frag_overlap'] = str(round(float(df[1][0].strip('%')),2)) +'%'
            stat['call_peak_number'] = df[1][1]
            stat['overlap_call_peak'] = str(round(float(df[1][2].strip('%')),2)) +'%'
            stat['percent_dup'] = str(round(float(df[1][3].strip('%')),2)) +'%'
        if i==9:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['nc_free_region'] = str(round(float(df[1][0].strip('%')),2)) +'%'
            stat['mono_nc_region'] = str(round(float(df[1][1].strip('%')),2)) +'%'
        if i==10:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['map_rate'] = df[1][1]
            stat['properly_reads'] = format(int(df[1][2]),",")
            stat['mit_rate'] = df[1][0]
        if i==11:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['rna_fastq_path'] = df[1][0]
            stat['ref_path_rna'] = df[1][1]
        if i==12:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['rna_saturation'] = df[1][0]
            stat['atac_saturation'] = df[1][1]
    plot_file = [
    path+'/report/div/barcode_rank.div',\
    path+'/report/div/cluster_chsize.div',\
    path+'/report/base64/6.base64',\
    path+'/report/base64/7.base64',\
    path+'/report/div/nUMI_chsize.div',\
    path+'/report/div/saturation.div',\
    path+'/report/div/saturation2.div',\
    path+'/report/div/anno_chsize.div',\
    path+'/report/div/barcode_rank_atac.div',\
    path+'/report/div/jaccard_rank.div',\
    path+'/report/base64/plot3_DropBeadsnum.base64',\
    path+'/report/base64/plot4_QC.base64',\
    path+'/report/base64/plot5_InterSize.base64',\
    path+'/report/base64/plot6_TSS.base64',\
    path+'/report/base64/plot7_Cluster_peak.base64',\
    path+'/report/base64/plot8_Cluster_depth.base64',\
    path+'/report/base64/plot9_Cluster_annotation.base64',\
    path+'/report/div/saturation_atac.div',\
    path+'/report/div/cluster_chsize-atac.div',\
    path+'/report/div/nUMI_chsize-atac.div',\
]
    plot_base64 = []
    plot_base64.append(open(path+'/report/div/barcode_rank.div',"r").read())
    plot_base64.append(open(path+'/report/div/cluster_chsize.div',"r").read())
    plot_base64.append(open(path+'/report/base64/6.base64',"r").read())
    plot_base64.append(open(path+'/report/base64/7.base64',"r").read())
    plot_base64.append(open(path+'/report/div/nUMI_chsize.div',"r").read())
    plot_base64.append(open(path+'/report/div/saturation.div',"r").read())
    plot_base64.append(open(path+'/report/div/saturation2.div',"r").read())
    plot_base64.append(open(path+'/report/div/barcode_rank_atac.div',"r").read())
    plot_base64.append(open(path+'/report/div/jaccard_rank.div',"r").read())
    plot_base64.append(open(path+'/report/base64/plot3_DropBeadsnum.base64',"r").read())
    plot_base64.append(open(path+'/report/base64/plot4_QC.base64',"r").read())
    plot_base64.append(open(path+'/report/base64/plot5_InterSize.base64',"r").read())
    plot_base64.append(open(path+'/report/base64/plot6_TSS.base64',"r").read())
    plot_base64.append(open(path+'/report/base64/plot7_Cluster_peak.base64',"r").read())
    plot_base64.append(open(path+'/report/base64/plot8_Cluster_depth.base64',"r").read())
    plot_base64.append(open(path+'/report/base64/plot9_Cluster_annotation.base64',"r").read())
    plot_base64.append(open(path+'/report/div/saturation_atac.div',"r").read())
    plot_base64.append(open(path+'/report/div/cluster_chsize-atac.div',"r").read())
    plot_base64.append(open(path+'/report/div/nUMI_chsize-atac.div',"r").read())

    if os.path.exists(path+'/report/div/anno_chsize.div'):
        plot_base64.append(open(path+'/report/div/anno_chsize.div',"r").read())

    plot_order = ['plot1','plot2','plot3','plot4','plot5','plot6','plot7','plot8','plot11','plot12','plot13','plot14','plot15','plot16','plot17','plot18','plot19','plot10','plot20','plot21']
    plot_dict = dict(zip(plot_order, plot_base64))
     
    data_tables_file = path+'/report/table/marker-table.txt'
    #data_tables = []
    #data_tables_dict = {}
    #for f_ in data_tables_file:
    table = str()
    if os.path.exists(data_tables_file):
        #if re.search('marker-table.txt',f_):
            table='''<table id=\"table_id_example\" class=\"table table-bordered table-striped\" style=\"Dosis;\">            <thead style=\"font-size:11px\"><tr>
                <th>gene</th>
                <th>cluster</th>
                <th>p_val_adj</th>
                <th>p_val</th>
                <th>avg_logFC</th>
                <th>pct.1</th>
                <th>pct.2</th>
                

                

            </tr>
        </thead>
            <tbody style=\"font-size:11px;\">
            '''+open(data_tables_file).read()+"</tbody></table>"
            #data_tables.append(table)               
        
    else:
        table.append('''
        <p style="font_family=DIN Next LT Pro;font_size=18px;font_weight=400">
        The table has not been generated because the data quality is too low.
        <p>
        ''')
    #print(data_tables)
    #table_order = ['table1','table2']
    #data_tables_dict = dict(zip(table_order, data_tables))
    import locale
    locale.setlocale(locale.LC_ALL, 'en_US')
    # for k,v in stat.items():
    #     if is_number(v):
    #         stat[k] =locale.format_string("%d", int(v), grouping=True)
    #     else:
    #         continue
    return stat, plot_dict, table
    
def write_param_to_template():
    stat, plot_dict, table = get_args_from_file()
    template = open(get_args()[1]).read()
    ID = get_args()[2]
    from string import Template
    path=get_args()[0]
    html=Template(template)
    if os.path.exists(path+'/report/div/anno_chsize.div'):
        report=html.safe_substitute(sample_info=ID, samplename=stat['samplename'],\
                    species=stat['species'], median_UMI_per_c=stat['median_UMI_per_c'],\
                    estm_Num_cell=stat['estm_Num_cell'],\
                    total_gene=stat['total_gene'],\
                    cluster_cell=stat['cluster_cell'],\
                       cDNA_num_frag=stat['cDNA_num_frag'],\
                       cDNA_frag_pass_QC=stat['cDNA_frag_pass_QC'],cDNA_frag_exact_bar=stat['cDNA_frag_exact_bar'],\
                       cDNA_frag_fail_bar=stat['cDNA_frag_fail_bar'],\
                       cDNA_frag_low_qual=stat['cDNA_frag_low_qual'],\
                       cDNA_Q30_c_bar=stat['cDNA_Q30_c_bar'],cDNA_Q30_s_bar=stat['cDNA_Q30_s_bar'],\
                       cDNA_Q30_UMI=stat['cDNA_Q30_UMI'],\
                       cDNA_Q30_r=stat['cDNA_Q30_r'],cDNA_adapter_filter=stat['cDNA_adapter_filter'],\
                       raw_r=stat['raw_r'],\
                       map_r=stat['map_r'],plus_strd=stat['plus_strd'],\
                       minus_strd=stat['minus_strd'],\
                       mito_ratio = stat['mito_ratio'],\
                       map_qual_corrt_r=stat['map_qual_corrt_r'],plot1=plot_dict['plot1'],\
                       plot2=plot_dict['plot2'],plot3=plot_dict['plot3'],\
                       plot4=plot_dict['plot4'],plot5=plot_dict['plot5'],plot6=plot_dict['plot6'],plot7=plot_dict['plot7'],plot8=plot_dict['plot8'],table = table,r_m_geno=stat['r_m_geno'],
                r_m_ex=stat['r_m_ex'],
                r_m_intro=stat['r_m_intro'],
                r_m_anti=stat['r_m_anti'],
                r_m_inter=stat['r_m_inter'],
                mean_r_per_c=stat['mean_r_per_c'],
                mean_UMI_per_c=stat['mean_UMI_per_c'],
                mean_genes_per_c=stat['mean_genes_per_c'],
                median_genes_per_c=stat['median_genes_per_c'],intron_boolean=stat['intron_boolean'],
                estimate_num_of_cell=stat['estimate_num_of_cell'],\
                    median_frag_per_cell=stat['median_frag_per_cell'], median_frac_peaks=stat['median_frac_peaks'],\
                    median_frac_tss=stat['median_frac_tss'],sample_id=stat['sample_id'],\
                    fastq_path=stat['fastq_path'],rna_fastq_path=stat['rna_fastq_path'], pipeline_version=stat['pipeline_version'],\
                    ref_path=stat['ref_path'],ref_path_rna=stat['ref_path_rna'],\
                    read_pairs=stat['read_pairs'],\
                    frac_valid_barcode=stat['frac_valid_barcode'],\
                    bead_thres=stat['bead_thres'],\
                    bead_number=stat['bead_number'],\
                    jaccard_thres=stat['jaccard_thres'],\
                    frac_frag_overlap=stat['frac_frag_overlap'],\
                    nc_free_region=stat['nc_free_region'],mono_nc_region=stat['mono_nc_region'],\
                    call_peak_number=stat['call_peak_number'],overlap_call_peak=stat['overlap_call_peak'],\
                    map_rate = stat['map_rate'],properly_reads = stat['properly_reads'],mit_rate = stat['mit_rate'],\
                    percent_dup=stat['percent_dup'],plot11=plot_dict['plot11'],\
                    plot12=plot_dict['plot12'],plot13=plot_dict['plot13'],\
                    plot14=plot_dict['plot14'],plot15=plot_dict['plot15'],\
                    plot16=plot_dict['plot16'],plot17=plot_dict['plot17'],\
                    plot18=plot_dict['plot18'],plot19=plot_dict['plot19'],\
              
                    plot10=plot_dict['plot10'],rna_saturation=stat['rna_saturation'],atac_saturation=stat['atac_saturation'],plot20=plot_dict['plot20'],plot21=plot_dict['plot21']
                
         
                
                )
           #plot6=plot_dict['plot6']
    else:
        report=html.safe_substitute(sample_info=ID, samplename=stat['samplename'],\
                    species=stat['species'], median_UMI_per_c=stat['median_UMI_per_c'],\
                    estm_Num_cell=stat['estm_Num_cell'],\
                    total_gene=stat['total_gene'],\
                    cluster_cell=stat['cluster_cell'],\
                       cDNA_num_frag=stat['cDNA_num_frag'],\
                       cDNA_frag_pass_QC=stat['cDNA_frag_pass_QC'],cDNA_frag_exact_bar=stat['cDNA_frag_exact_bar'],\
                       cDNA_frag_fail_bar=stat['cDNA_frag_fail_bar'],\
                       cDNA_frag_low_qual=stat['cDNA_frag_low_qual'],\
                       cDNA_Q30_c_bar=stat['cDNA_Q30_c_bar'],cDNA_Q30_s_bar=stat['cDNA_Q30_s_bar'],\
                       cDNA_Q30_UMI=stat['cDNA_Q30_UMI'],cDNA_adapter_filter=stat['cDNA_adapter_filter'],\
                       cDNA_Q30_r=stat['cDNA_Q30_r'],\
                       raw_r=stat['raw_r'],\
                       map_r=stat['map_r'],plus_strd=stat['plus_strd'],\
                       minus_strd=stat['minus_strd'],\
                       mito_ratio = stat['mito_ratio'],\
                       map_qual_corrt_r=stat['map_qual_corrt_r'],plot1=plot_dict['plot1'],\
                       plot2=plot_dict['plot2'],plot3=plot_dict['plot3'],\
                       plot4=plot_dict['plot4'],plot5=plot_dict['plot5'],table = table,r_m_geno=stat['r_m_geno'],
                       plot6=plot_dict['plot6'],plot7=plot_dict['plot7'],plot8="There is no such species reference for annnotation.",
                r_m_ex=stat['r_m_ex'],
                r_m_intro=stat['r_m_intro'],
                r_m_anti=stat['r_m_anti'],
                r_m_inter=stat['r_m_inter'],
                mean_r_per_c=stat['mean_r_per_c'],
                mean_UMI_per_c=stat['mean_UMI_per_c'],
                mean_genes_per_c=stat['mean_genes_per_c'],
                median_genes_per_c=stat['median_genes_per_c'],intron_boolean=stat['intron_boolean'],
                
                estimate_num_of_cell=stat['estimate_num_of_cell'],\
                    median_frag_per_cell=stat['median_frag_per_cell'], median_frac_peaks=stat['median_frac_peaks'],\
                    median_frac_tss=stat['median_frac_tss'],sample_id=stat['sample_id'],\
                    fastq_path=stat['fastq_path'],rna_fastq_path=stat['rna_fastq_path'], pipeline_version=stat['pipeline_version'],\
                    ref_path=stat['ref_path'],ref_path_rna=stat['ref_path_rna'],\
                    read_pairs=stat['read_pairs'],\
                    frac_valid_barcode=stat['frac_valid_barcode'],\
                    bead_thres=stat['bead_thres'],\
                    bead_number=stat['bead_number'],\
                    jaccard_thres=stat['jaccard_thres'],\
                    frac_frag_overlap=stat['frac_frag_overlap'],\
                    nc_free_region=stat['nc_free_region'],mono_nc_region=stat['mono_nc_region'],\
                    call_peak_number=stat['call_peak_number'],overlap_call_peak=stat['overlap_call_peak'],\
                    map_rate = stat['map_rate'],properly_reads = stat['properly_reads'],mit_rate = stat['mit_rate'],\
                    percent_dup=stat['percent_dup'],plot11=plot_dict['plot11'],\
                    plot12=plot_dict['plot12'],plot13=plot_dict['plot13'],\
                    plot14=plot_dict['plot14'],plot15=plot_dict['plot15'],\
                    plot16=plot_dict['plot16'],plot17=plot_dict['plot17'],\
                    plot18=plot_dict['plot18'],plot19="There is no such species reference for annnotation.",\
                    #plot9=plot_dict['plot9'],
                    plot10=plot_dict['plot10'],rna_saturation=stat['rna_saturation'],atac_saturation=stat['atac_saturation'],plot20=plot_dict['plot20'],plot21=plot_dict['plot21']
                
                
                
                
                
                
                
                
                
                
                
                
                )
           #plot6=plot_dict['plot6']        

    return report
    
if __name__ == '__main__':
    outpath=get_args()[0]
    ID = get_args()[2]
    x = datetime.datetime.now()
    dd = x.strftime("%Y%m%d")
    get_args_from_file()
    report=write_param_to_template()
    fw = open(outpath+'/report/'+ID+'_'+dd+'_C4-scCAT_report.html','w')
    fw.write(report)
