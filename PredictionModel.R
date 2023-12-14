### Load packages
library(GenomicFeatures)
library(ChIPseeker)
library(randomForest)
library(caret)
library(class)
library(rtracklayer)

### Validation of the model using test Beko data
Beko_test_data <- readRDS("x_rf.rds")
rf_classifier <- readRDS("rf_Model.rds")
accuracy <- function(x) {sum(diag(x)/(sum(rowSums(x)))*100)}
nor <- function(x) {(x-min(x))/(max(x)-min(x))}

### Data structure
head(Beko_test_data)

### Preidction of dynamic and static sites
prediction_for_table <- predict(rf_classifier,Beko_test_data)
tab <- table(observed=Beko_test_data[,1],predicted=prediction_for_table)
tab
accuracy(tab)

### Generate own dataset
### Load MSPC output (including cooridnates and MSPC's p-value)
cDat <- read.table("./ConsensusPeaks.bed", header = T)
cDat_GR <- GRanges(cDat[1:3])

### Include -log10(FDR)
cDat_GR$pval <- cDat$X.1xlog10.p.value.
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
cDat_GR$Distance <- NA
cDat_GR[paste(cDat_GR) %in% paste0(peakAnno_cDat$seqnames,":", peakAnno_cDat$start,"-",peakAnno_cDat$end)]$Distance <- peakAnno_cDat$annotation
cDat_GR <- cDat_GR[!is.na(cDat_GR$Distance)]

cDat_DF <- as.data.frame(cDat_GR)

rownames(cDat_DF) <- paste(cDat_GR)

cDat_DF$Distance <- sapply(strsplit(cDat_DF$Distance , "Intron "), "[", 1)
cDat_DF$Distance [is.na(cDat_DF$Distance )] <- "Intron"

cDat_DF$Distance <- sapply(strsplit(cDat_DF$Distance , "Exon "), "[", 1)
cDat_DF$Distance [is.na(cDat_DF$Distance )] <- "Exon"
cDat_DF$Distance [cDat_DF$Distance  == ""] <- "Exon"

### Remove currently not included features
#cDat_DF <- cDat_DF[cDat_DF$Distance != "Downstream (<=300bp)",]
cDat_DF <- cDat_DF[cDat_DF$Distance %in% Beko_test_data$Distance,]

### Convert into factors
cDat_DF$Distance <- as.factor(cDat_DF$Distance)

### Temporary fix: Include missing Factors from Test_data
test <- rbind(cDat_DF[,7:8], Beko_test_data[,-1])
prediction_for_table <- predict(rf_classifier,test)

### Remove added factors
prediction_for_table <- prediction_for_table[1:nrow(cDat_DF)]

### Prediction outcome
table(prediction_for_table)

### Export .bed files
Dynamic_Sites <- GRanges(names(prediction_for_table[prediction_for_table == "Dynamic"]))
Static_Sites <- GRanges(names(prediction_for_table[prediction_for_table == "Static"]))

export.bed(Dynamic_Sites, "Predicted_Dynamic_Sites.bed")
export.bed(Static_Sites, "Predicted_Static_Sites.bed")
