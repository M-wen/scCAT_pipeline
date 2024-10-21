library(Seurat)
library(data.table)

parser = argparse::ArgumentParser(description="")
parser$add_argument('--matrix', help='RAW MEX matrix dir')
parser$add_argument('--hex',help='16 2 ATCG file')
parser$add_argument('--tran', help='ATAC beads2cell file')
parser$add_argument('--output',help='Output dir, default: current dir')
args = parser$parse_args()



tran <- fread(args$tran, header=F)
hex <- fread(args$hex, header=F)
data <- Read10X(args$matrix,gene.column=1)


tmp <- as.data.frame(colSums(data))
tmp$BARCODE<- rownames(tmp)
colnames(tmp) <- c("UMI","BARCODE")


dd <- merge(tmp, hex, by.x="BARCODE", by.y="V1")
all <- merge(dd, tran, by.x="V2", by.y="V1", all=T)
colnames(all)[4] <- "Beads"

all[which(is.na(all$Beads)),"Beads"] <- 'noise'
all[which(!(all$Beads %in% "noise")),"Beads"] <- "true"
alll<- all[which(!(is.na(all$UMI)) ),]

alls <- alll[order(alll$UMI, decreasing=T),c("UMI","Beads")]

alls$barcodes <- c(1:dim(alls)[1])

write.csv(alls[,c("barcodes","UMI","Beads")], paste0(args$output,"/2.cutoff.csv"),quote=F, row.names=F)
