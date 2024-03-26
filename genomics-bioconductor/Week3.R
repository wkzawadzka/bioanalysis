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
