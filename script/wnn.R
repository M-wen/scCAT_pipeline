library(harmony)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(data.table)
library(GenomicRanges)
library(stringr)
library(scMCA)
library(scHCL)
library(rhdf5)
library(svglite)
library(tidyr)
library(SoupX)
library(DropletUtils)
library(dplyr)
library(cowplot)
library(DoubletFinder)
library(RColorBrewer)
#library(hdf5r)
set.seed(1234)

### Get the parameters
parser = argparse::ArgumentParser(description="merging wnn, QC and BC function")
parser$add_argument('-R','--rna', help='input rna FilterMatrix file')
parser$add_argument('-A','--atac', help='input atac Peak file')
parser$add_argument('-SN','--samplename', help='sample name default:None')
parser$add_argument('-F','--frag', help='input *fragments.tsv.gz')
parser$add_argument('-SP','--species', help='hg38 or mm10')
parser$add_argument('-RM','--RawMatrix', help='input rna RawMatrix file')
parser$add_argument('-BT','--BarcodeTran', help='input barcodeTranslate_16.txt file')
parser$add_argument('-D','--dim',help='dim usage')
parser$add_argument('-MP','--mtgenepercentage',help='filter cells with mtgenes percentage')
parser$add_argument('-MF','--minfeatures',help='filter cells with minimum nfeatures')
parser$add_argument('-PC','--pc',help='pc usage')
parser$add_argument('-RES','--res',help='resolution usage')
parser$add_argument('-K','--knn',help='defines k for the k-nearest neighbor algorithm')
parser$add_argument('-MD','--maxdim',help='max dimension to keep from UMAP procedure') 
parser$add_argument('-O','--outdir', help='outputdir default:None')
parser$add_argument('-IS','--ifSoupx', help='whether use SoupxMatrix to analysis')
parser$add_argument('-MT','--chrmt', help='type of mt')
parser$add_argument('-GTF',"--gtf", help='gtf file')

args = parser$parse_args()
species.usage <- if(!is.null(args$species)) args$species else "mm10"

# if (species.usage == "mouse"||species.usage == "Mouse"||species.usage == "mm10"||species.usage == "Mouse_nucleus" || species.usage == "Mouse_V2.3"){
#   gtf <- "/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/database/Mouse/genes.gtf"
# }else if (species.usage == "human"||species.usage == "Human"||species.usage == "hg19"||species.usage == "hg38" ||species.usage == "Human_nucleus" || species.usage == "hg19_premRNA" || 
#          species.usage == "Human_2020A_mkgtf" || species.usage == "Human_2020A" ||species.usage == "Human_93" || species.usage == "Human_V2.3" ||species.usage == "Human_optimizedv1"){
#   gtf <- "/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/database/Human/genes.gtf" 
# }
gtf <- if(!is.null(args$gtf)) args$gtf else "/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/database/Mouse/genes.gtf"


###########################
# overlap and save h5file #
###########################
# rna
rna_counts <- Read10X(args$rna,gene.column = 1)
# atac
mtx_path <- paste(args$atac, "/matrix.mtx", sep = '')
feature_path <- paste(args$atac, "/peak.bed", sep = '')
barcode_path <- paste(args$atac, "/barcodes.tsv", sep = '')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
atac_counts <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

# intersect
cells<- intersect(colnames(atac_counts),colnames(rna_counts))
rna_counts_cells <- rna_counts[,which(colnames(rna_counts) %in% cells)]
atac_counts_cells <- atac_counts[,which(colnames(atac_counts) %in% cells)]


# save.h5
saveh5file <- function(rna_counts_cells,atac_counts_cells){
  cell_name <- rna_counts_cells@Dimnames[[2]]
  multi.data <- rbind(rna_counts_cells,atac_counts_cells)
  h5createFile(paste(args$outdir,"/overlap_matrix.h5",sep = ''))
  # Saving matrix information.
  h5createGroup(paste(args$outdir,"/overlap_matrix.h5",sep = ''),"matrix")
  h5write(multi.data@Dimnames[[2]] , paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/barcodes")
  h5write(multi.data@x, paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/data")
  
  h5createGroup(paste(args$outdir,"/overlap_matrix.h5",sep = ''),"matrix/features")
  key <- c('genome','interval')
  h5write(key, paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/features/_all_tag_keys")
  Genes <- rep('Gene Expression', length(rna_counts_cells@Dimnames[[1]]))
  Peaks <- rep("Peaks", length(atac_counts_cells@Dimnames[[1]]))
  Features <- c(Genes,Peaks)
  h5write(Features,paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/features/feature_type")
  Genome <- rep(args$species, length(multi.data@Dimnames[[1]]))
  h5write(Genome,paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/features/genome")
  h5write(multi.data@Dimnames[[1]],paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/features/id")
  
  h5write(multi.data@Dimnames[[1]],paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/features/name")
  h5write(multi.data@i, paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/indices") # already zero-indexed.
  h5write(multi.data@p, paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/indptr")
  h5write(dim(multi.data), paste(args$outdir,"/overlap_matrix.h5",sep = ''), "matrix/shape")
  h5closeAll()
}
unlink(paste(args$outdir,"/overlap_matrix.h5",sep = ''),recursive=T,force=T)
saveh5file(rna_counts_cells,atac_counts_cells)

print('save.h5 file has already!')


#########
# SoupX #
#########
run_soupx <- function(toc,tod,rho=NULL) {
  tod <- tod[rownames(toc),]
  all <- toc
  all <- CreateSeuratObject(all)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(all)
  all <- ScaleData(all, features = all.genes)
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
  all <- FindNeighbors(all, dims = 1:30)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:30)
  
  matx <- all@meta.data
  
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
  sc = setContaminationFraction(sc, 0.2)
  out = adjustCounts(sc)
  saveRDS(sc,"sc.rds")
  DropletUtils:::write10xCounts(paste(args$outdir,"/../RNA/03.Matrix/SoupxMatrix",sep=""), out,version="3")
}

toc <- Read10X(args$rna,gene.column = 1)
bac <- Read10X(args$RawMatrix,gene.column = 1)
valid <- read.table(args$BarcodeTran,header = F)
bac <- bac[rownames(toc),]
bac <- bac[,which(!(colnames(bac) %in% valid$V1))]
tod <- cbind(bac,toc)
unlink(paste(args$outdir,"/../RNA/03.Matrix/SoupxMatrix",sep=""),recursive=T,force=T)
run_soupx(toc,tod)

print('soupx has already!')

########################
# User choosing matrix #
########################
if(args$ifSoupx %in% c('t','T','true','True')){
  object <- Read10X(paste(args$outdir,"/../RNA/03.Matrix/SoupxMatrix",sep=""),gene.column=1)
  object <- CreateSeuratObject(object)
}else{
  object <- CreateSeuratObject(rna_counts_cells)
}


#################
# QC_analysis.R #
#################
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(as.numeric(doublets.percentage)*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  #data<-subset(data,subset=doublet_info=="Singlet")
  data
}
nfeature_plot <- function(EC){
  extra_nfeature <- data.frame(FetchData(object = EC, vars = 'nFeature_RNA'))
  extra_nfeature$group <- "nFeature_RNA"
  p <- ggplot(extra_nfeature, aes(x = group, y = nFeature_RNA, fill=group)) + 
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
    ylab(expression("nFeature_RNA"))+xlab("")
  ylims_feature <- extra_nfeature %>%
    group_by(extra_nfeature$group) %>%
    summarise(Q1 = quantile(extra_nfeature$nFeature_RNA, 1/4,na.rm=T), Q3 = quantile(extra_nfeature$nFeature_RNA, 3/4,na.rm=T)) %>%
    ungroup() %>%
    #get lowest Q1 and highest Q3
    summarise(lowQ1 = 0, highQ3 = max(Q3)*4)
  p <- p + coord_cartesian(ylim = as.numeric(ylims_feature))
  p
}
ncount_plot <- function(EC){
  extra_ncount <- data.frame(FetchData(object = EC, vars = 'nCount_RNA'))
  extra_ncount$group <- "nCount_RNA"
  p <- ggplot(extra_ncount, aes(x = group, y = nCount_RNA, fill=group)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F,outlier.fill = "white",outlier.colour = "white") + 
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
    ylab(expression("nCount_RNA"))+xlab("")
  ylims_feature <- extra_ncount %>%
    group_by(extra_ncount$group) %>%
    summarise(Q1 = quantile(extra_ncount$nCount_RNA, 1/4,na.rm=T), Q3 = quantile(extra_ncount$nCount_RNA, 3/4,na.rm=T)) %>%
    ungroup() %>%
    #get lowest Q1 and highest Q3
    summarise(lowQ1 = 0, highQ3 = max(Q3)*4)
  p <- p + coord_cartesian(ylim = as.numeric(ylims_feature))
  p
}
percentMt_plot <- function(EC){
  extra_mt <- data.frame(FetchData(object = EC, vars = 'percent.mt'))
  extra_mt$group <- "percent.mt"
  p <- ggplot(extra_mt, aes(x = group, y = percent.mt, fill=group)) + 
    geom_violin(trim=TRUE,color="white",show.legend = F,outlier.fill = "white",outlier.colour = "white") + 
    geom_boxplot(width=0.06,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0, 
                 outlier.stroke = 0)+ 
    scale_fill_manual(values = "#56B4E9")+
    theme_cowplot()+ 
    theme(axis.text.x=element_blank(), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(family="Times",size=14,face="plain"), 
          axis.title.y=element_text(family="Times",size = 18,face="plain"))+
    guides(fill="none")+
    ylab(expression("percent.mt"))+xlab("")
  ylims_feature <- extra_mt %>%
    group_by(extra_mt$group) %>%
    summarise(Q1 = quantile(extra_mt$percent.mt, 1/4,na.rm=T), Q3 = quantile(extra_mt$percent.mt, 3/4,na.rm=T)) %>%
    ungroup() %>%
    #get lowest Q1 and highest Q3
    summarise(lowQ1 = 0, highQ3 = max(Q3)*20)
  p <- p + coord_cartesian(ylim = as.numeric(ylims_feature))
  p
}


dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 20)
doublets.percentage <- as.numeric(if(!is.null(args$percentage)) args$percentage else 0.05)
mtgene_path <- if(!is.null(args$mtgenes)) args$mtgenes else "auto"
mtegne_filter <- as.numeric(if(!is.null(args$mtgenepercentage)) args$mtgenepercentage else 10)
minfeatures <- as.numeric(if(!is.null(args$minfeatures)) args$minfeatures else 200)

### Creat Seurat object
grange.counts <- StringToGRanges(rownames(atac_counts_cells), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts_cells <-  atac_counts_cells[as.vector(grange.use), ]

gtf <- rtracklayer::import(gtf)
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

if(species.usage %in% c("mouse","Mouse","mm10","human","Human","hg19","hg38","Human_nucleus","hg19_premRNA","Mouse_nucleus","Human_2020A_mkgtf","Human_2020A","Human_93","Human_V2.3","Human_optimizedv1","Mouse_V2.3")){
  gene.coords$gene_biotype <- gene.coords$gene_type
}
  
frag.file <- args$frag
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts_cells,
  sep = c(":", "-"),
  fragments = frag.file,
  min.cells = 10,
  annotation = gene.coords
)

object[["ATAC"]] <- chrom_assay
DefaultAssay(object) <- "ATAC"
object <- TSSEnrichment(object)
total_fragments <- CountFragments(frag.file)
rownames(total_fragments)<- total_fragments$CB
object@meta.data$fragments <- total_fragments[colnames(object), "frequency_count"]
object <- FRiP(
  object = object,
  assay = 'ATAC',
  total.fragments = 'fragments'
)
object@meta.data$loguniqueFrag <- log10(object@meta.data$fragments)

######
# QC #
######
DefaultAssay(object) <- "RNA"
object <- NormalizeData(object)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
object <- ScaleData(object)

if(dim(object)[2] <=50){
  object <- RunPCA(object,npcs = (dim(object)[2]-1))
}else{
  object <- RunPCA(object)
}
object <- RunUMAP(object, dims = 1:dim.usage,reduction.name = "umap.rna")

mtgene <- length(grep ("^mt-", rownames(object[["RNA"]]),value = T))
print(mtgene)
MTgene <- length(grep ("^MT-", rownames(object[["RNA"]]),value = T))
print(MTgene)

#### Plot raw and filter QC Vlnplot
png(paste(args$outdir,"/plot/",args$samplename,"_raw_QCplot.png",sep=""))
if(mtgene_path != "auto"){
  mt_gene_table <- read.table(mtgene_path,sep="\t")
  mtgene <- as.character(mt_gene_table[,1])
  object[["percent.mt"]] <- PercentageFeatureSet(object, features = mtgene)
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p3 <- percentMt_plot(object)
  p <- p1|p2|p3
}else if(mtgene>0){
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p3 <- percentMt_plot(object)
  p <- p1|p2|p3
}else if(MTgene>0){
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p3 <- percentMt_plot(object)
  p <- p1|p2|p3
}else{
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p <- p1|p2
}
print(p)
dev.off()
ggsave(paste(args$outdir,"/plot/",args$samplename,"_raw_QCplot.svg",sep=""), p , width = 7, height = 4)


### Filter cells with nfeatures/percent.mt
png(paste(args$outdir,"/plot/",args$samplename,"_filter_QCplot.png",sep=""))
objectmeta <- object@meta.data[order(-object@meta.data$nFeature_RNA),]
n95 <- as.numeric(as.integer(nrow(objectmeta) * 0.05))
n95features <- as.numeric(objectmeta[n95,"nFeature_RNA"])
if(mtgene_path != "auto"){
  object <- subset(object, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features & percent.mt < mtegne_filter)
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p3 <- percentMt_plot(object)
  p <- p1|p2|p3
}else if(mtgene>0){
  object <- subset(object, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features & percent.mt < mtegne_filter)
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p3 <- percentMt_plot(object)
  p <- p1|p2|p3
}else if(MTgene>0){
  object <- subset(object, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features & percent.mt < mtegne_filter)
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p3 <- percentMt_plot(object)
  p <- p1|p2|p3
}else{
  object <- subset(object, subset = nFeature_RNA > minfeatures & nFeature_RNA < n95features)
  p1 <- nfeature_plot(object)
  p2 <- ncount_plot(object)
  p <- p1|p2
}
print(p)
dev.off()
ggsave(paste(args$outdir,"/plot/",args$samplename,"_filter_QCplot.svg",sep=""), p , width = 7, height = 4)

#Find doublets
if(dim(object)[2] >50){
  object <- Find_doublet(object)
  write.table(object@meta.data,paste0(args$outdir,"/",args$samplename,"_doublets_info.txt"),sep="\t",quote=FALSE)
  object <- subset(object,subset=doublet_info=="Singlet")
}
object@meta.data$split = args$samplename

print('QC plotting has already!')


############
# analysis #
############

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(object) <- "ATAC"
object <- RunTFIDF(object)
object <- FindTopFeatures(object, min.cutoff = 'q0')
object <- RunSVD(object)
object <- RunUMAP(object, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

object <- FindMultiModalNeighbors(object, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
object <- RunUMAP(object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(object, reduction = "umap.rna", label = TRUE,pt.size = 0.8, label.size = 4, repel = TRUE) + ggtitle("RNA")+NoLegend()
p2 <- DimPlot(object, reduction = "umap.atac", label = TRUE,pt.size = 0.8, label.size = 4, repel = TRUE) + ggtitle("ATAC")+NoLegend()
p3 <- DimPlot(object, reduction = "wnn.umap", label = TRUE,pt.size = 0.8, label.size = 4, repel = TRUE) + ggtitle("WNN")
ggsave(paste(args$outdir,"/plot/wnn_cluster.png",sep = ''),p1+p2+p3,width=26,height=8)
ggsave(paste(args$outdir,"/plot/wnn_cluster.svg",sep = ''),p1+p2+p3,width=26,height=8)

print('wnn cluster plotting has already!')
saveRDS(object,paste0(args$outdir,"/WNN_cluster.RDS"))


## get clutering data for html plotting
cluster_ID=as.data.frame(Idents(object = object))
cluster_cor= as.data.frame(Embeddings(object = object,reduction = "umap.atac"))
coor=cbind(cluster_ID,cluster_cor,object[['nCount_ATAC']],object[['nFeature_ATAC']])
colnames(coor) = c("Cluster","UMAP_1","UMAP_2","nCount_ATAC","nFeature_ATAC")
coorOrder = coor[order(coor$Cluster),]

temp <- coorOrder
names <- rownames(temp)
rownames(temp) <- NULL
dataTemp <- cbind(names,temp)

rm(temp,names,coorOrder,coor,cluster_cor,cluster_ID)

cluster_stat <- as.data.frame(table(dataTemp$Cluster))
colnames(cluster_stat) <- c("Cluster","cellNum")
cluster_cell <- dplyr::left_join(dataTemp,cluster_stat,by="Cluster")

length = nrow(cluster_stat)
write.csv(cluster_stat, file=paste(args$outdir,"/cluster_cell_atac.stat",sep=""),quote=FALSE)


# atac annotation
Annotation(object) <- gene.coords
gene.activities <- GeneActivity(object)

if (species.usage %in% c("mouse","Mouse","mm10","human","Human","hg19","hg38","Human_nucleus","hg19_premRNA","Mouse_nucleus","Human_2020A_mkgtf","Human_2020A","Human_93","Human_V2.3","Human_optimizedv1","Mouse_V2.3")){
  if (species.usage == "mouse"||species.usage == "Mouse"||species.usage == "mm10"||species.usage == "Mouse_nucleus" || species.usage == "Mouse_V2.3"){
    object.combined_cell_type <- scMCA(scdata = gene.activities, numbers_plot = 3)
    out=as.data.frame(unlist(object.combined_cell_type$scMCA))
    out$`unlist(object.combined_cell_type$scMCA)`=as.character(out$`unlist(object.combined_cell_type$scMCA)`)
  }
  else if (species.usage == "human"||species.usage == "Human"||species.usage == "hg19"||species.usage == "hg38" ||species.usage == "Human_nucleus" || species.usage == "hg19_premRNA" || species.usage == "Human_2020A_mkgtf" || species.usage == "Human_2020A" ||species.usage == "Human_93" || species.usage == "Human_V2.3" ||species.usage == "Human_optimizedv1"){
    object.combined_cell_type <- scHCL(scdata = gene.activities, numbers_plot = 3)
    out=as.data.frame(unlist(object.combined_cell_type$scHCL))
    out$`unlist(object.combined_cell_type$scHCL)`=as.character(out$`unlist(object.combined_cell_type$scHCL)`)
  }
  object@meta.data$cell_type_atac=out[match(rownames(object@meta.data),rownames(out)),1]
  out_meta=object@meta.data

  table_list=list()
  for(i in 0:(length(unique(object$seurat_clusters))-1)){
    a=i+1
    table_list[[a]] <- tryCatch(
      {
        sub=subset(out_meta,out_meta$seurat_clusters==i)
        tab=as.data.frame(table(sub$cell_type_atac))
        tab_order=tab[order(tab[,2],decreasing = T),]
        tab_order$Cluster=i
        #table_list[[a]]=tab_order[1,]
        tab_order[1,]
      },error=function(e){table_list[[a]]='NA'}
    )
  }
  sum=do.call(rbind,table_list)
  use=out_meta[,c("seurat_clusters","cell_type_atac")]
  use$seurat_clusters=as.character(use$seurat_clusters)
  use$ID=rownames(out_meta)
  colnames(sum)=c("predicated.cell.type","Freq","seurat_clusters")
  res=merge(use,sum,by="seurat_clusters")
  object@meta.data$predicated.cell.type.atac=res[match(rownames(out_meta),res$ID),4]
}else{
  print('This refcode has not auto annotation')
}

## create cluster.csv
if(species.usage %in% c("mouse","Mouse","mm10","human","Human","hg19","hg38","Human_nucleus","hg19_premRNA","Mouse_nucleus","Human_2020A_mkgtf","Human_2020A","Human_93","Human_V2.3","Human_optimizedv1","Mouse_V2.3")){
  cluster_anno <- as.data.frame(table(object@meta.data$cell_type_atac,object@meta.data$seurat_clusters))
  sorted_cluster_anno <- cluster_anno[order(cluster_anno$Var2,-cluster_anno$Freq),]
  finl_cluster_anno <- sorted_cluster_anno[!duplicated(sorted_cluster_anno$Var2),]
  finl_cluster_anno <- finl_cluster_anno[,c(1,2)]

  colnames(finl_cluster_anno) <- c("CellType","Cluster")
  finl_cluster_anno$CellType <- as.character(finl_cluster_anno$CellType)

  celltype_cell <- dplyr::left_join(cluster_cell,finl_cluster_anno,by="Cluster")
  temp <- as.data.frame(table(celltype_cell$CellType))
  colnames(temp) <- c("CellType","CellTypeNum")
  celltype_cell <- dplyr::left_join(celltype_cell,temp,by="CellType")
  rm(temp)
  celltype_cell_merge <- unite(celltype_cell, "Cluster", Cluster, cellNum, sep = " CellsNum: ")
  celltype_cell_merge <- unite(celltype_cell_merge, "Predicted cell type", CellType, CellTypeNum, sep = ": ")
  rownames(celltype_cell_merge) <- celltype_cell_merge[,1]
  celltype_cell_merge <- celltype_cell_merge[,-1]

  write.csv(celltype_cell_merge, file=paste(args$outdir,"/report/ATAC/6.cluster.csv",sep=""),quote=FALSE)
  
  # plot
  p1 <- DimPlot(object, reduction = "umap.rna", group.by = "predicated.cell.type.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")+NoLegend()
  p2 <- DimPlot(object, reduction = "umap.atac", group.by = "predicated.cell.type.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")+NoLegend()
  p3 <- DimPlot(object, reduction = "wnn.umap", group.by = "predicated.cell.type.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  ggsave(paste(args$outdir,"/plot/wnn_annotation_atac.png",sep = ''),p1+p2+p3,width=26,height=8)
  ggsave(paste(args$outdir,"/plot/wnn_annotation_atac.svg",sep = ''),p1+p2+p3,width=26,height=8)
  
}else{
  object@meta.data$cell_type_atac <- 'NULL'
  finl_cluster_anno <- as.data.frame(table(object@meta.data$cell_type_atac,object@meta.data$seurat_clusters))
  finl_cluster_anno <- finl_cluster_anno[,c(1,2)]
  colnames(finl_cluster_anno) <- c("CellType","Cluster")
  finl_cluster_anno$CellType <- as.character(finl_cluster_anno$CellType)
  
  celltype_cell <- dplyr::left_join(cluster_cell,finl_cluster_anno,by="Cluster")
  temp <- as.data.frame(table(celltype_cell$CellType))
  colnames(temp) <- c("CellType","CellTypeNum")
  celltype_cell <- dplyr::left_join(celltype_cell,temp,by="CellType")
  rm(temp)
  celltype_cell_merge <- unite(celltype_cell, "Cluster", Cluster, cellNum, sep = " CellsNum: ")
  celltype_cell_merge <- unite(celltype_cell_merge, "Predicted cell type", CellType, CellTypeNum, sep = ": ")
  rownames(celltype_cell_merge) <- celltype_cell_merge[,1]
  celltype_cell_merge <- celltype_cell_merge[,-1]
  
  write.csv(celltype_cell_merge, file=paste(args$outdir,"/report/ATAC/6.cluster.csv",sep=""),quote=FALSE)
  
  # plot
  p1 <- DimPlot(object, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")+NoLegend()
  p2 <- DimPlot(object, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")+NoLegend()
  p3 <- DimPlot(object, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  ggsave(paste(args$outdir,"/plot/wnn_annotation_atac.png",sep = ''),p1+p2+p3,width=26,height=8)
  ggsave(paste(args$outdir,"/plot/wnn_annotation_atac.svg",sep = ''),p1+p2+p3,width=26,height=8)
}

print('atac annotation has already!')



# rna annotation
DefaultAssay(object) <- "RNA"
rna_data <- GetAssayData(object = object, slot = "data")

if (species.usage %in% c("mouse","Mouse","mm10","human","Human","hg19","hg38","Human_nucleus","hg19_premRNA","Mouse_nucleus","Human_2020A_mkgtf","Human_2020A","Human_93","Human_V2.3","Human_optimizedv1","Mouse_V2.3")){
  if (species.usage == "mouse"||species.usage == "Mouse"||species.usage == "mm10"||species.usage == "Mouse_nucleus" || species.usage == "Mouse_V2.3"){
    object.combined_cell_type <- scMCA(scdata = rna_data, numbers_plot = 3)
    out=as.data.frame(unlist(object.combined_cell_type$scMCA))
    out$`unlist(object.combined_cell_type$scMCA)`=as.character(out$`unlist(object.combined_cell_type$scMCA)`)
  }
  else if (species.usage == "human"||species.usage == "Human"||species.usage == "hg19"||species.usage == "hg38" ||species.usage == "Human_nucleus" || species.usage == "hg19_premRNA" || species.usage == "Human_2020A_mkgtf" || species.usage == "Human_2020A" ||species.usage == "Human_93" || species.usage == "Human_V2.3" ||species.usage == "Human_optimizedv1"){
    object.combined_cell_type <- scHCL(scdata = rna_data, numbers_plot = 3)
    out=as.data.frame(unlist(object.combined_cell_type$scHCL))
    out$`unlist(object.combined_cell_type$scHCL)`=as.character(out$`unlist(object.combined_cell_type$scHCL)`)
  }

  object@meta.data$cell_type_rna=out[match(rownames(object@meta.data),rownames(out)),1]
  out_meta=object@meta.data

  table_list=list()
  for(i in 0:(length(unique(object$seurat_clusters))-1)){
  a=i+1
  table_list[[a]] <- tryCatch(
    {
      sub=subset(out_meta,out_meta$seurat_clusters==i)
      tab=as.data.frame(table(sub$cell_type_rna))
      tab_order=tab[order(tab[,2],decreasing = T),]
      tab_order$Cluster=i
      tab_order[1,]
    },error=function(e){table_list[[a]]='NA'}
  )
  sum=do.call(rbind,table_list)
  use=out_meta[,c("seurat_clusters","cell_type_rna")]
  use$seurat_clusters=as.character(use$seurat_clusters)
  use$ID=rownames(out_meta)
  colnames(sum)=c("predicated.cell.type","Freq","seurat_clusters")
  res=merge(use,sum,by="seurat_clusters")
  object@meta.data$predicated.cell.type.rna=res[match(rownames(out_meta),res$ID),4]
  }
  
  print('rna annotation has already!')

  # plot
  p1 <- DimPlot(object, reduction = "umap.rna", group.by = "predicated.cell.type.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")+NoLegend()
  p2 <- DimPlot(object, reduction = "umap.atac", group.by = "predicated.cell.type.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")+NoLegend()
  p3 <- DimPlot(object, reduction = "wnn.umap", group.by = "predicated.cell.type.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  ggsave(paste(args$outdir,"/plot/wnn_annotation_rna.png",sep = ''),p1+p2+p3,width=26,height=8)
  ggsave(paste(args$outdir,"/plot/wnn_annotation_rna.svg",sep = ''),p1+p2+p3,width=26,height=8)
  
}else{
  p1 <- DimPlot(object, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")+NoLegend()
  p2 <- DimPlot(object, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")+NoLegend()
  p3 <- DimPlot(object, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  ggsave(paste(args$outdir,"/plot/wnn_annotation_rna.png",sep = ''),p1+p2+p3,width=26,height=8)
  ggsave(paste(args$outdir,"/plot/wnn_annotation_rna.svg",sep = ''),p1+p2+p3,width=26,height=8)
  
  print('This refcode has not auto annotation')
}

##-----------------------------
## create cluster_cell_rna.stat
cluster_ID=as.data.frame(Idents(object = object))
cluster_cor= as.data.frame(Embeddings(object = object,reduction = "umap.rna"))
coor=cbind(cluster_ID,cluster_cor,object[['nCount_RNA']],object[['nFeature_RNA']])
colnames(coor) = c("Cluster","UMAP_1","UMAP_2","nUMI","nGene")
coorOrder = coor[order(coor$Cluster),]
temp <- coorOrder
names <- rownames(temp)
rownames(temp) <- NULL
dataTemp <- cbind(names,temp)

rm(temp)
rm(names)
rm(coorOrder)
rm(coor)
rm(cluster_cor)
rm(cluster_ID)

cluster_stat <- as.data.frame(table(dataTemp$Cluster))
colnames(cluster_stat) <- c("Cluster","cellNum")
cluster_cell <- dplyr::left_join(dataTemp,cluster_stat,by="Cluster")
write.csv(cluster_stat, file=paste(args$out,"/cluster_cell_rna.stat",sep=""),quote=FALSE)

##-----------------------
## create cluster_rna.csv
if(species.usage %in% c("mouse","Mouse","mm10","human","Human","hg19","hg38","Human_nucleus","hg19_premRNA","Mouse_nucleus","Human_2020A_mkgtf","Human_2020A","Human_93","Human_V2.3","Human_optimizedv1","Mouse_V2.3")){
  cluster_anno <- as.data.frame(table(object@meta.data$cell_type_rna,object@meta.data$seurat_clusters))
  sorted_cluster_anno <- cluster_anno[order(cluster_anno$Var2,-cluster_anno$Freq),]
  finl_cluster_anno <- sorted_cluster_anno[!duplicated(sorted_cluster_anno$Var2),]
  finl_cluster_anno <- finl_cluster_anno[,c(1,2)]
  
  colnames(finl_cluster_anno) <- c("CellType","Cluster")
  finl_cluster_anno$CellType <- as.character(finl_cluster_anno$CellType)
  
  celltype_cell <- dplyr::left_join(cluster_cell,finl_cluster_anno,by="Cluster")
  temp <- as.data.frame(table(celltype_cell$CellType))
  colnames(temp) <- c("CellType","CellTypeNum")
  celltype_cell <- dplyr::left_join(celltype_cell,temp,by="CellType")
  rm(temp)
  celltype_cell_merge <- unite(celltype_cell, "Cluster", Cluster, cellNum, sep = " CellsNum: ")
  celltype_cell_merge <- unite(celltype_cell_merge, "Predicted cell type", CellType, CellTypeNum, sep = ": ")
  rownames(celltype_cell_merge) <- celltype_cell_merge[,1]
  celltype_cell_merge <- celltype_cell_merge[,-1]
  
  write.csv(celltype_cell_merge, file=paste(args$outdir,"/cluster_rna.csv",sep=""),quote=FALSE)
}else{
  object@meta.data$cell_type_rna <- 'NULL'
  finl_cluster_anno <- as.data.frame(table(object@meta.data$cell_type_rna,object@meta.data$seurat_clusters))
  finl_cluster_anno <- finl_cluster_anno[,c(1,2)]
  colnames(finl_cluster_anno) <- c("CellType","Cluster")
  finl_cluster_anno$CellType <- as.character(finl_cluster_anno$CellType)
  
  celltype_cell <- dplyr::left_join(cluster_cell,finl_cluster_anno,by="Cluster")
  temp <- as.data.frame(table(celltype_cell$CellType))
  colnames(temp) <- c("CellType","CellTypeNum")
  celltype_cell <- dplyr::left_join(celltype_cell,temp,by="CellType")
  rm(temp)
  celltype_cell_merge <- unite(celltype_cell, "Cluster", Cluster, cellNum, sep = " CellsNum: ")
  celltype_cell_merge <- unite(celltype_cell_merge, "Predicted cell type", CellType, CellTypeNum, sep = ": ")
  rownames(celltype_cell_merge) <- celltype_cell_merge[,1]
  celltype_cell_merge <- celltype_cell_merge[,-1]
  
  write.csv(celltype_cell_merge, file=paste(args$outdir,"/cluster_rna.csv",sep=""),quote=FALSE)
}

saveRDS(object,paste0(args$outdir,"/WNN_annotation.RDS"))


# FindAllMarker
allmarkers.rna <- FindAllMarkers(object, assay = 'RNA')
write.csv(as.data.frame(allmarkers.rna[,c(6,5,1,2,3,4)]),file= paste0(args$outdir,"/",args$samplename,"_rna_marker.csv"),quote=FALSE)

cat(paste("Number of cells used for clustering,", length(colnames(object)), "\n",sep=""),
    file=paste(args$outdir,"/cell_report_2.csv",sep=""))


# statistic of the mitochondria mapping rate
fragment <- as.data.frame(fread(frag.file,header=F))
MT <- fragment[grep(paste0("^",args$chrmt),fragment$V1),]
m1=data.frame(qc="Mitochondria reads ratio",num=paste0(100*as.numeric(sprintf("%0.4f",sum(MT$V5)/sum(fragment$V5))),"%"),stringsAsFactors = FALSE)
write.table(m1,paste0(args$outdir,"/report/ATAC/3.mapping.csv"),sep=":",quote = FALSE,row.names = FALSE,col.names = FALSE)


print('WNN.R part has already!')
