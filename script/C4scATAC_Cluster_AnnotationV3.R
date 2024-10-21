### Get the parameters
parser = argparse::ArgumentParser(description="Script to clustering and celltype annotation scATAC data")
parser$add_argument('-I','--peak', help='input narrowpeak file')
parser$add_argument('-F','--fragment', help='input fragment files')
parser$add_argument('-C','--cell', help='input cell list')
parser$add_argument('-G','--promoter', help='input promoter bed file')
parser$add_argument('-Q','--qc', help='input qc file')
parser$add_argument('-O','--out', help='out directory')
parser$add_argument('-MT','--chrmt', help='chrmt')
args = parser$parse_args()


library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(data.table)
library(Matrix)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(stringr)
library(scMCA)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

peak_set <- as.data.frame(fread(args$peak,sep="\t",header = F))
peak_set <- peak_set[,c(1,2,3)]
colnames(peak_set) <- c("chr", "start", "end")
grfile <- makeGRangesFromDataFrame(peak_set)

cells <- read.table(args$cell,header=F)
frag <- CreateFragmentObject(path=args$fragment,cells=cells$V1)
counts <- FeatureMatrix(fragments = frag,features = grfile, cells=cells$V1,sep=c(":","-"))

## write peak-cell matrix to Matrix Market Formats
writeMM(counts, file=paste0(args$out,"/01.out/Peak/matrix.mtx"))
write.table(colnames(counts),file=paste0(args$out,"/01.out/Peak/barcodes.tsv"),quote=F,row.names=F,col.names=F)
write.table(rownames(counts),file=paste0(args$out,"/01.out/Peak/peak.bed"),quote=F,row.names=F,col.names=F)


## QC 
counts <- counts[grep(paste0("^",args$chrmt),rownames(counts),invert=T),]
raw_meta=read.table(args$qc,header=T)
rownames(raw_meta)=raw_meta[,1]
raw_meta=raw_meta[,-1]
my_meta=raw_meta
colnames(my_meta)[2]="uniqueNuclearFrags" #0604 modified
#my_meta$FRIP=colSums(counts)/2/my_meta$uniqueNuclearFrags ## 2022-05-24 modified V3 matrix count fragment number not divided by 2
my_meta$FRIP=colSums(counts)/my_meta$uniqueNuclearFrags
my_meta$log10_uniqueFrags=log10(my_meta$uniqueNuclearFrags)
my_meta$Sample="scATAC"

QC_plot=list()
QC_plot[[1]] <- ggplot(my_meta, aes(x = Sample, y = log10_uniqueFrags, fill=Sample)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F,outlier.fill = "white",outlier.colour = "white") + 
    geom_boxplot(width=0.06,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0, 
    outlier.stroke = 0)+ 
    scale_fill_manual(values = "#999999")+
    theme_cowplot()+ 
    theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family="Times",size=14,face="plain"), 
        axis.title.y=element_text(family="Times",size = 18,face="plain"))+
    guides(fill="none")+
    ylab(expression(Log[10]*paste("(","Fragments",")",sep = "")))+xlab("") 

QC_plot[[2]] <- ggplot(my_meta, aes(x = Sample, y = tssProportion, fill=Sample)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F) + 
    geom_boxplot(width=0.06,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0,
    outlier.stroke = 0)+ 
    scale_fill_manual(values = "#E69F00")+
    theme_cowplot()+ 
    theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family="Times",size=14,face="plain"), 
        axis.title.y=element_text(family="Times",size = 18,face="plain"))+
    guides(fill="none")+
    ylab("TSS Proportion")+xlab("") 

QC_plot[[3]] <- ggplot(my_meta, aes(x = Sample, y = FRIP, fill=Sample)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F) + 
    geom_boxplot(width=0.06,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0,
    outlier.stroke = 0)+ 
    scale_fill_manual(values = "#56B4E9")+
    theme_cowplot()+ 
    theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family="Times",size=14,face="plain"), 
        axis.title.y=element_text(family="Times",size = 18,face="plain"))+
    guides(fill= "none")+
    ylab("FRIP")+xlab("") 

ggsave(paste0(args$out,"/../Joint/report/ATAC/plot4_QC.png"), do.call(plot_grid,c(QC_plot, ncol = 3, align = "hv")), width = 7, height = 4)
ggsave(paste0(args$out,"/../Joint/report/ATAC/plot4_QC.svg"), do.call(plot_grid,c(QC_plot, ncol = 3, align = "hv")), width = 7, height = 4)

###output report file

qc5=data.frame(qc="Fraction of fragments overlapping TSS",num=paste0(100*round(sum(my_meta$tssProportion*my_meta$uniqueNuclearFrags)/sum(my_meta$uniqueNuclearFrags),5),"%"),stringsAsFactors = FALSE)
qc5[2,1]="Called peak number"
qc5[2,2]=prettyNum(nrow(counts),big.mark = ",")
qc5[3,1]="Fraction of fragments overlapping called peaks"
qc5[3,2]=paste0(100*round(sum(my_meta$FRIP*my_meta$uniqueNuclearFrags)/sum(my_meta$uniqueNuclearFrags),5),"%")
qc5[4,1]="Percent duplicates"
qc5[4,2]=paste0(100*round(1-(sum(my_meta$uniqueNuclearFrags)/sum(my_meta$totalFrags)),5),"%")
write.table(qc5,paste0(args$out,"/../Joint/report/ATAC/5.library.QC.csv"),sep = ":",quote = FALSE,row.names = FALSE,col.names = FALSE)

qc1=data.frame(qc="Estimated number of cells",num=prettyNum(ncol(counts),big.mark = ","),stringsAsFactors = FALSE)
qc1[2,1]="Median fragments per cell"
qc1[2,2]=prettyNum(median(my_meta$uniqueNuclearFrags),big.mark = ",")
qc1[3,1]="Median fraction of fragments overlapping peaks"
qc1[3,2]=paste0(100*round(median(my_meta$FRIP),5),"%")
qc1[4,1]="Median fraction of fragments overlapping TSSs"
qc1[4,2]=paste0(100*round(median(my_meta$tssProportion),5),"%")
write.table(qc1,paste0(args$out,"/../Joint/report/ATAC/1.cell_report.csv"),sep = ":",quote = FALSE,row.names = FALSE,col.names = FALSE)

###Number of beads per droplet

name=as.data.frame(rownames(my_meta))
colnames(name)="DropBarcode"
name$Num=unlist(lapply(strsplit(as.character(name$DropBarcode),"N0"),"[",2))
table=as.data.frame(table(name$Num))
colnames(table)=c("Num","Count")

png(paste0(args$out,"/../Joint/report/ATAC/plot3_DropBeadsnum.png"),width = 458,height = 377)
p=ggplot(data = table, mapping = aes(x = factor(Num), y = Count, fill = Num)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_brewer(palette = 'Set2',labels = paste(levels(table$Num)," ",table$Count))+
  theme_bw()+
  xlab("Number of beads per droplet")+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  ggtitle(paste("Total cell number",nrow(my_meta)))
print(p)
dev.off()

svg(paste0(args$out,"/../Joint/report/ATAC/plot3_DropBeadsnum.svg"),width = 5,height = 4)
print(p)
dev.off()

### clustering
#assay <- CreateChromatinAssay(counts, fragments = frag,sep = c(":","-"))
#scATAC <- CreateSeuratObject(assay, assay = 'peaks',meta.data = my_meta)
scATAC <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'C4scATAC',
  min.cells = 10,
  meta.data = my_meta
)

scATAC <- subset(scATAC, subset = uniqueNuclearFrags > 3 & tssProportion > 0.1)
scATAC <- RunTFIDF(scATAC)
scATAC <- FindTopFeatures(scATAC, min.cutoff = 'q0')
scATAC <- RunSVD(
  object = scATAC,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
scATAC <- RunUMAP(object = scATAC, reduction = 'lsi', dims = 1:30)
scATAC <- FindNeighbors(object = scATAC, reduction = 'lsi', dims = 1:30)
scATAC <- FindClusters(object = scATAC, verbose = FALSE)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
a=getPalette(length(unique(scATAC@meta.data$seurat_clusters)))

c1=DimPlot(object = scATAC, label = TRUE) + NoLegend()+
  scale_color_manual(values = a)

png(paste0(args$out,"/../Joint/report/ATAC/plot7_Cluster_peak.png"),width = 524,height = 488)
print(c1)
dev.off()

svg(paste0(args$out,"/../Joint/report/ATAC/plot7_Cluster_peak.svg"),width = 5,height = 5)
print(c1)
dev.off()

c2=FeaturePlot(scATAC,features = "log10_uniqueFrags")+
  scale_color_viridis(direction = -1,option = "D",name="log10 Frags")+
  ggtitle("")

png(paste0(args$out,"/../Joint/report/ATAC/plot8_Cluster_depth.png"),width = 524,height = 488)
print(c2)
dev.off()

svg(paste0(args$out,"/../Joint/report/ATAC/plot8_Cluster_depth.svg"),width = 6,height = 5)
print(c2)
dev.off()

###Cell annotation

## write promoter-cell matrix to Matrix Market Formats
#tss <- as.data.frame(fread("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/USER/liuyang9/pipe/C4scATAC/database/refdata_wxy_Macaca_fascicularis/Promoter_4k_Macaca_fascicularis.bed",sep="\t",header = F))
tss=read.table(args$promoter,header = F)
## Jul.20,2022 remove the chromesome of promoter bed files not present in fragment files
fragment <- as.data.frame(fread(args$fragment,header=F))
tss <- tss[which(tss$V1 %in% unique(fragment$V1)),]
## Jul.25,2022 statistic of the mitochondria mapping rate 
MT <- fragment[grep(paste0("^",args$chrmt),fragment$V1),]
m1=data.frame(qc="Mitochondria reads ratio",num=paste0(100*as.numeric(sprintf("%0.4f",sum(MT$V5)/sum(fragment$V5))),"%"),stringsAsFactors = FALSE)
write.table(m1,paste0(args$out,"/../Joint/report/ATAC/3.mapping.csv"),sep = ":",quote = FALSE,row.names = FALSE,col.names = FALSE)
rm(MT)
rm(m1)

###
tss$Name=unlist(lapply(strsplit(as.character(tss$V4),":"),"[",1))
tss$Name=make.unique(tss$Name)
tss$Name=str_to_title(tss$Name)
#rownames(tss) <- make.unique(tss$Name)
tmp <- tss[,c(1,2,3)]
colnames(tmp) <- c("chr", "start", "end")
grpromoter <- makeGRangesFromDataFrame(tmp)
promoter <- FeatureMatrix(fragments = frag,features = grpromoter, cells=cells$V1,sep=c(":","-"))
writeMM(promoter, file=paste0(args$out,"/01.out/Promoter/matrix.mtx"))
write.table(colnames(promoter),file=paste0(args$out,"/01.out/Promoter/barcodes.tsv"),quote=F,row.names=F,col.names=F)
write.table(rownames(promoter),file=paste0(args$out,"/01.out/Promoter/promoter.bed"),quote=F,row.names=F,col.names=F)
rm(tmp)

rownames(promoter) <- make.unique(tss$Name)
scATAC_result <- scMCA(scdata = promoter, numbers_plot = 3)
if (length(scATAC_result$scMCA) != 0){
  out=as.data.frame(unlist(scATAC_result$scMCA))
  out$`unlist(scATAC_result$scMCA)`=as.character(out$`unlist(scATAC_result$scMCA)`)

  scATAC@meta.data$cell_type=out[match(rownames(scATAC@meta.data),rownames(out)),1]

  out_meta=scATAC@meta.data

  table_list=list()

  for(i in 0:(length(unique(scATAC$seurat_clusters))-1)){
    a=i+1
    sub=subset(out_meta,out_meta$seurat_clusters==i)
    tab=as.data.frame(table(sub$cell_type))
    tab_order=tab[order(tab[,2],decreasing = T),]
    tab_order$Cluster=i
    table_list[[a]]=tab_order[1,]
  }

  sum=do.call(rbind,table_list)
  use=out_meta[,c(15,16)]
  use$seurat_clusters=as.character(use$seurat_clusters)
  use$ID=rownames(out_meta)
  colnames(sum)=c("predicated.cell.type","Freq","seurat_clusters")
  res=merge(use,sum,by="seurat_clusters")
  scATAC@meta.data$predicated.cell.type=res[match(rownames(out_meta),res$ID),4]

  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  a=getPalette(length(unique(scATAC@meta.data$predicated.cell.type)))

  c3=DimPlot(object = scATAC, label = FALSE, group.by = "predicated.cell.type") +
    scale_color_manual(values = a)

  png(paste0(args$out,"/../Joint/report/ATAC/plot9_Cluster_annotation.png"),width = 990,height = 479)
  print(c3)
  dev.off()

  svg(paste0(args$out,"/../Joint/report/ATAC/plot9_Cluster_annotation.svg"),width = 11,height = 5)
  print(c3)
  dev.off()

  table_cell.type=as.data.frame(table(as.character(scATAC@meta.data$predicated.cell.type)))
  colnames(table_cell.type)=c("predicated.cell.type","number")
  table_cell.type$ratio=table_cell.type$number/colSums(table_cell.type[2])
  order_table_cell.type=table_cell.type[order(table_cell.type$number,decreasing = T),]
  write.table(order_table_cell.type,paste0(args$out,"/../Joint/report/ATAC/Table2.celltypecount.csv"),sep = ",",quote = FALSE,row.names = FALSE)

} else{
  table_cell.type=data.frame(predicated.cell.type="None",number=dim(scATAC)[2],ratio=1)
  write.table(table_cell.type,paste0(args$out,"/../Joint/report/ATAC/Table2.celltypecount.csv"),sep = ",",quote = FALSE,row.names = FALSE)
}
saveRDS(scATAC,paste0(args$out,"/../Joint/report/ATAC/saved_clustering.rds"))

###Find differentially accessible peaks between clusters
# da_peaks=FindAllMarkers(scATAC,
#                         min.pct = 0.2,
#                         test.use = 'LR')
# sub_da_peak=subset(da_peaks,da_peaks$p_val_adj<0.01 & da_peaks$avg_log2FC > log(exp(1)))
# sub_da_peak$chr=unlist(lapply(strsplit(sub_da_peak$gene,":"),"[",1))
# sub_da_peak$tag=unlist(lapply(strsplit(sub_da_peak$gene,":"),"[",2))
# sub_da_peak$start=unlist(lapply(strsplit(sub_da_peak$tag,"-"),"[",1))
# sub_da_peak$end=unlist(lapply(strsplit(sub_da_peak$tag,"-"),"[",2))
# sub_da_peak=sub_da_peak[,-9]
# 
# use=sub_da_peak[,c(8,9,10,7)]
# colnames(use)[4]="peak ID"
# tss=tss[,c(1:3,5)]
# colnames(tss)=c("chr","start","end","gene symbol")
# GRange_TSS=makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)
# GRange_Peak=makeGRangesFromDataFrame(use, keep.extra.columns = TRUE)
# distance=distanceToNearest(GRange_Peak,GRange_TSS)
# out_Peak=as.data.frame(GRange_Peak[queryHits(distance)])[,c(6,4)]
# out_TSS_GRange <- GRange_TSS[subjectHits(distance)]
# 
# out_TSS=as.data.frame(out_TSS_GRange)[6]
# summary=cbind(out_Peak,out_TSS)
# summary$distance=distance@elementMetadata$distance
# sub_da_peak=sub_da_peak[,c(1:7)]
# colnames(sub_da_peak)[7]="peak.ID"
# result_da_peak=unique(merge(summary,sub_da_peak,by="peak.ID"))
# result_da_peak=result_da_peak[order(result_da_peak$cluster),]
# write.table(result_da_peak,paste0(args$out,"/../Joint/report/ATAC/Table1.dapeak.csv"),sep = ",",quote = FALSE,row.names = FALSE)

