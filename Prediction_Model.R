### Load packages
library(GenomicFeatures)
library(ChIPseeker)
library(randomForest)
library(caret)
library(class)
library(rtracklayer)

### Validation of the model using test data set
test_data_set <- readRDS("test_data_rf.rds")
rf_classifier <- readRDS("rf_Model.rds")
accuracy <- function(x) {sum(diag(x)/(sum(rowSums(x)))*100)}
nor <- function(x) {(x-min(x))/(max(x)-min(x))}

### Data structure
head(test_data_set)

### Preidction of dynamic and static sites
prediction_for_table <- predict(rf_classifier,test_data_set[,-1])
tab <- table(observed=test_data_set[,1],predicted=prediction_for_table)
tab
accuracy(tab)

### Generate own dataset
### Load MSPC output (including cooridnates and MSPC's p-value)
cDat <- read.table("./New_system/RBPJ_Ref.bed", header = F)
cDat_GR <- GRanges(seqnames = cDat$V1, IRanges(start = cDat$V2, end = cDat$V3))

### Include -log10(FDR)
cDat_GR$pval <- cDat$V5
### min-max normalization of MSPC's p-vlaues
xnorm <- nor(cDat_GR$pval)

cDat_GR$MSPC_FDR <- xnorm


### Load txdb of Genome of Interest based on gtf
#txdb <- makeTxDbFromGFF("/home/tobiasf/Genomes/hg19/hg19_genes_iGenomes.gtf", format = "gtf", organism = "Homo sapiens")

### Save and load txdb
#saveDb(txdb, "txdb_hg19")
txdb <- loadDb("txdb_hg19")
### Annotation of genomic features using annotatePeak function
peakAnno_cDat <- annotatePeak(cDat_GR, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno_cDat <- as.data.frame(peakAnno_cDat)

### Remove missing annotation
cDat_GR$Feature <- NA
cDat_GR[paste(cDat_GR) %in% paste0(peakAnno_cDat$seqnames,":", peakAnno_cDat$start,"-",peakAnno_cDat$end)]$Feature <- peakAnno_cDat$annotation
cDat_GR <- cDat_GR[!is.na(cDat_GR$Feature)]

cDat_DF <- as.data.frame(cDat_GR)
cDat_DF <- cDat_DF[,c(7,8)]
rownames(cDat_DF) <- paste(cDat_GR)



### Reduce features and combine Introns and Exons
cDat_DF$Feature <- sapply(strsplit(cDat_DF$Feature , "Intron "), "[", 1)
cDat_DF$Feature [is.na(cDat_DF$Feature )] <- "Exon"
cDat_DF$Feature [cDat_DF$Feature  == ""] <- "Exon"

cDat_DF$Feature <- sapply(strsplit(cDat_DF$Feature , "Exon "), "[", 1)
cDat_DF$Feature [is.na(cDat_DF$Feature )] <- "Exon"
cDat_DF$Feature [cDat_DF$Feature  == ""] <- "Exon"

cDat_DF$Feature <- sapply(strsplit(cDat_DF$Feature , "Downstream "), "[", 1)
cDat_DF$Feature[cDat_DF$Feature  == ""] <- "Downstream"
cDat_DF$Feature[is.na(cDat_DF$Feature)] <- "Downstream"

cDat_DF$Feature <- sapply(strsplit(cDat_DF$Feature , "Promoter "), "[", 1)
cDat_DF$Feature[cDat_DF$Feature  == ""] <- "Promoter"
cDat_DF$Feature[is.na(cDat_DF$Feature)] <- "Promoter"



### Include binding motifs using FIMO output
RBPJ_Stat <- read.delim("New_system/RBPJ_stat.tsv")
RBPJ_StatGR <- GRanges(RBPJ_Stat$sequence_name)

RBPJ_Dyn <- read.delim("New_system/RBPJ_dyn.tsv")
RBPJ_DynGR <- GRanges(RBPJ_Dyn$sequence_name)

RBPJ_SP <- read.delim("New_system/RBPJ_stat.tsv")
RBPJ_SPGR <- GRanges(RBPJ_SP$sequence_name)



cDat_DF$Motif_RBPJDyn <- ifelse(GRanges(rownames(cDat_DF)) %over% RBPJ_DynGR, "True","False")
cDat_DF$Motif_RBPStat <- ifelse(GRanges(rownames(cDat_DF)) %over% RBPJ_StatGR, "True","False")
cDat_DF$Motif_SP <- ifelse(GRanges(rownames(cDat_DF)) %over% RBPJ_SPGR, "True","False")

### Convert into factors
cDat_DF$Feature <- as.factor(cDat_DF$Feature)
cDat_DF$Motif_RBPJDyn <- as.factor(cDat_DF$Motif_RBPJDyn)
cDat_DF$Motif_RBPStat <- as.factor(cDat_DF$Motif_RBPStat)
cDat_DF$Motif_SP <- as.factor(cDat_DF$Motif_SP)


### Prediction of the generated data frame
prediction_for_table <- predict(rf_classifier,cDat_DF)


### Export .bed files
Dynamic_Sites <- GRanges(names(prediction_for_table[prediction_for_table == "Dynamic"]))
Static_Sites <- GRanges(names(prediction_for_table[prediction_for_table == "Static"]))

export.bed(Dynamic_Sites, "Predicted_Dynamic_Sites.bed")
export.bed(Static_Sites, "Predicted_Static_Sites.bed")

