## Week 4 
##### Illumina
# Taught by: Kasper Daniel Hansen, PhD
# -------- 450K DNA methylation array ---------
# DNA methylation is a chemical modification of
# the C base that in humans only occurs in a 
# CpG context. (C followed by G in the human genome)

# Popular platform for starting DNA methylation
# is the 450K microarray - comprehesive and cheap
# we have ∼28million CpG's

library(minfi)
library(GEOquery)

# study: wheather or not DNA methylation changes
# are associated with acute mania
# # of indiv. hospitalized for acute mania -> obtained serum
# from them -> profiled the serum on the 450K array
accesion_number = "GSE68777"
getGEOSuppFiles(accesion_number)

untar(paste0(accesion_number,"/", accesion_number, "_RAW.tar"), exdir=paste0(accesion_number, "/raw"))
list.files(paste0(accesion_number, "/raw")) # idat.gz"

# two files for each sample: green & red

# decompress
idatFiles = list.files(paste0(accesion_number, "/raw"), full.names = TRUE)
sapply(idatFiles, gunzip, overwrite=TRUE)

# read expression data
RGSet <- read.metharray.exp(paste0(accesion_number, "/raw"))
RGSet # annotation: ilmn12.hg19

assay(RGSet)[1:10]
# [1]  234 6197 2220 1326 2514 2749  974  689 3687  151
max(assay(RGSet)) # [1] 21772

# get Phenotype info
mat <- getGEO(accesion_number)
mat[[1]]
pD.all <- pData(mat[[1]])
head(pD.all)
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)

# clean up names
names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ","", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)
sampleNames(RGSet) <- sub(".*_5", "5", sampleNames(RGSet))
rownames(pD) <- pD$title

# merge
pD <- pD[sampleNames(RGSet), ] # reorder
head(pD)
pDd <- as(pD, "DataFrame")
pData(RGSet) <- pDd
#showClass("RGChannelSet")
head(RGSet)

# preprocess
GRSet <- preprocessQuantile(RGSet)
?preprocessQuantile
GRSet # normalized data & maps data to genome
# assign each probe to given location in genome that
# particular CpG is located

head(getIslandStatus(GRSet))
getBeta(GRSet)[1:3,1:3]

# we are interested in clusters of CpG's that change
# in the same diirection


### Count-based RNA-seq analysis
# ferature: gene or an exon
# we know: # num of reads matching a particular feature
library(DESeq2)
library(edgeR)
library(airway)

data(airway)
assay(airway, "counts")[1:3,1:3]
airway$dex #Levels: trt untrt

airway$dex <- relevel(airway$dex, "untrt")
airway$dex # Levels: untrt trt

# Why? THe reference level for factor is the first level (this case treatment)
# goal: find differences from treated -> untreated (a bit unnatural)

granges(airway)
dge <- DGEList(counts=assay(airway, "counts"), 
               group= airway$dex)
dge$samples <- merge(dge$samples,
                     as.data.frame(colData(airway)),
                     by=0)
head(dge$genes)
dge$genes <- data.frame(name = names(rowRanges(airway)),
                        stringsAsFactors = FALSE)
head(dge$genes)

# norm factors
dge <- calcNormFactors(dge)

# 
dge <- estimateGLMCommonDisp(dge)
dge <- estimateGLMTagwiseDisp(dge)
dge

# 
design <- model.matrix(~dge$samples$group)
fit <-glmFit(dge, design)

# 
lrt <- glmLRT(fit, coef = 2)
topTags(lrt) # most diff expressed


### DESeq2setup
library(DESeq2)
dds <- DESeqDataSet(airway, design = ~ dex)
# fit
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(red$padj),]
res

