# course: Bioconductor for Genomic Data Science
# Taught by: Kasper Daniel Hansen, PhD
### Week 3

library(Biobase)
BiocManager::install("ALL")
library(ALL)
BiocManager::install("hgu95av2.db")
library(hgu95av2.db)

data(ALL)
ALL

experimentData(ALL)
exprs(ALL)[1:4, 1:4]
head(pData(ALL)$sex)


## gene annotation mapping
ids = featureNames(ALL)[1:5]
ids
as.list(hgu95av2ENTREZID[ids])


phenoData(ALL)

# summarized experiment
BiocManager::install("airway")
library(airway)
library(GenomicRanges)
data(airway)
colData(airway)
head(rownames(airway))


assayNames(airway)
assay(airway, "counts")[1:4, 1:4]
length(rowRanges(airway))

rowRanges(airway)
# GRanges object with 17 ranges
# -> has 17 exons and is located on X chromosome


### GEO query
library(GEOquery)

eList <- getGEO("GSE11675")
class(eList)
length(eList)
names(eList)
eData <- eList[[1]]
eData
names(pData(eData))


eList_raw = getGEOSuppFiles("GSE11675")
eList_raw
untar("C:/Users/weraz/Documents/GSE11675/GSE11675_RAW.tar")


### biomaRt
library(biomaRt)
head(listMarts())
mart <- useMart("ensembl")
mart
head(listDatasets(mart))
ensembl <- useDataset("hsapiens_gene_ensembl", mart)



### R S4 Classes
## classes & methods
df <- data.frame(y = rnorm(10), x= rnorm(10))
lm.object <- lm(y~x, data=df)
lm.object
names(lm.object)
class(lm.object)

xx = list(a = letters[1:3], b=rnorm(4))
class(xx) = "lm"
xx


library(ALL)
data(ALL)
ALL # Expression Set
class(ALL)
isS4(ALL) # true
class?ExpressionSet


new("ExpressionSet") # constructor
getClass("ExpressionSet") # definition of a class
ALL@annotation #  "hgu95av2"
annotation(ALL)


# R S4 Methods
library(GRanges)

as.data.frame
base::as.data.frame
showMethods("as.data.frame")
getMethod("as.data.frame", signature(x = "GenomicRanges"))

method?"as.data.frame,DataFrame"

showMethods("findOverlaps")
getMethod("findOverlaps", signature(query = "Ranges", subject="Ranges"))
           



# test
ALL[, 5]
exprs(ALL[,5])
mean(exprs(ALL[,5]))


# annotate ALL
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
mart <- useMart(host='https://feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)

features <- featureNames(ALL) 
# getBM Retrieves information from the BioMart database
ann <- getBM(attributes = c("ensembl_gene_id", "affy_hg_u95av2"),
      filters    = "affy_hg_u95av2",
      values     = features, 
      mart       = ensembl)


# How many probesets (features) are annotated with more than one Ensembl gene id?
ann_names = ann[,2]
ann_names
sum(table(ann_names)>1)


# How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chromosomes 1 to 22).
table(ann[,1])
filters = listFilters(ensembl)
filters[1:5,]
chrom <- c(1:22)
ann_chrom <- getBM(attributes = c("ensembl_gene_id", "affy_hg_u95av2"),
             filters    = c("affy_hg_u95av2", "chromosome_name"),
             values     = list(features, chrom), 
             mart       = ensembl)
sum(table(table(ann_chrom[,2])))

# Question: What is the mean value of the Methylation channel across the features for sample “5723646052_R04C01”?
BiocManager::install("minfiData")
library(minfiData)
# MsetEx dataset
data(MsetEx)

MsetEx
assayNames(MsetEx)
mean(assay(MsetEx, "Meth")[, 2])

# Question: Access the processed data from NCBI GEO Accession number GSE788. What is the mean expression level of sample GSM9024?
library(GEOquery)
geo = getGEO("GSM9024")
mean(geo@dataTable@table$VALUE)


# Question: What is the average of the average length across the samples in the expriment?
library(airway)
data(airway)
mean(airway$avgLength)

# 7. Question 7 We are using the airway dataset from the airway package. The features in this dataset are Ensembl genes.
# What is the number of Ensembl genes which have a count of 1 read or more in sample SRR1039512?
# SRR1039512
colnames(airway)
airway[, 3]
as <- assay(airway[, 3])
sum(as>=1)

# 
chrom = c(1:22)
chrom
chrom_names = paste0("chr", chrom)

# How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22) as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene) #you will need to make sure that the transcript representation does not contain introns
exons

exons_autosome <- keepSeqlevels(exons, autosome, pruning.mode="coarse")

ncbi <- mapSeqlevels(seqlevels(exons_autosome), "NCBI")
exons_autosome_ncbi <- renameSeqlevels(exons_autosome, ncbi)

dim(subsetByOverlaps(airway, exons_autosome_ncbi))[1]
