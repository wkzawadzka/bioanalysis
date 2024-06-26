# course: Bioconductor for Genomic Data Science
# Taught by: Kasper Daniel Hansen, PhD
library(oligo)
library(GEOquery)
# Affymetrix: very cheap, very high quality, very short articles
# the probes in the array are typically on the order of 25 bases long
# thats not very long, it means the probes are not very specific, to compensate that
# on most affymetrix arrays  if you are measuring one specific on a species, or you're measuring a slip, you do this
# using multiple probes that all measure the same tack
# a probe set is a group of probes that all measure the same tack.
# [quote from Kasper Daniel Hansen, PhD]

# For technical reasons it's quite expensive to design an isometrics chip, but then once it's designed,
# it's cheap to mass produce

# These multiple probes are grouped into something called probe sets
# The use of (PM,MM) pairs and multiple pairs for a target transcript allows both absolute and comparative analysis
# and compensates for variations and noises in the complex background
# : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3902802/
# Affymetrix uses one-color method 

# Affymetrix microarray - get Data
accession = "GSE38792"
getGEOSuppFiles(accession)
list.files(accession)
untar(paste0(accession,"/", accession, "_RAW.tar"), exdir=paste0(accession, "/CEL"))
list.files(paste0(accession, "/CEL")


celfiles <- list.files(paste0(accession, "/CEL"), full = TRUE)
celfiles
rawData <- read.celfiles(celfiles)

# show
# experiment: comparing samples from some control patients
# to patiens with sleep apnea
rawData 
# more than 1 million probes! (1102500 features)

# pd.hugene.1.0.st.v1 -< based on random priming instead of OligoDT priming
# they have a lot of features for each gene
# human gene

# The features, the probes that you are using to base your
# RNA transcipt are spaced along the entire length od the transcript

getClass("GeneFeatureSet") # [package "oligoClasses"]

# we have many probes that measure a single gene
# and this all is taken care of in gene feature set


### let's look into the data
exprs(rawData)[1:4, 1:3]
# pretty large integers visible (200-10000)
# raw intensity measurement from the scanner

# so micro array scanner is typicall a 16 bit scanner
# which means when you scan a probe you get a number between
# zero and 2^16 = 65,536, to confirm:

max(exprs(rawData))
# [1] 65534

# research shows, this is not very nice range to work on.
# so we tranform it
# we take log2 and we get a number between zero and 16 :)

# pData
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)), "OSA", "Control")
pData(rawData)


# boxplot
boxplot(rawData, "all")
# y axis on the log scale
# so differene of 1 or 2 is quite massive difference!

# samples 5-7: low intensities and low spread
# hypothesis: nothing was hybridized for those samples


# Normalization RMA
normData <- rma(rawData)
# RMA:
# (1) background correction
# (2) quantile normalization 
# (3) take all diff probes that make up the same
# gene and output single number for that
normData
exprs(normData)[1:4, 1:3] # now we see log ranges
max(exprs(normData)) #[1] 14.37372


dim(rawData)
dim(normData)
# Huge reduction in features 1102500 -> 33297
# bacause of (3)rd step of RMA

# Isometrics identifiers
featureNames(normData)[1:10]
# [1] "7892501" "7892502" "7892503" "7892504"
# We need to go from these numbers to gene names now.

boxplot(normData, "all") # ready for analysis :)
# these boxplots are not identical because we run quantile normalization
# on the probe level and then we summarize from probe level into the gene level
# same if (1) summ -> then (2) norm


### diff expr
library(limma)

# model with contrasts
design2 <- model.matrix(~ normData$group -1)
colnames(design2) <- c("Control", "OSA")
# Control is the reference level.
contrasts <- makeContrasts("OSA-Control", levels=design2)
fit2 <- lmFit(normData, design2)
fit2C <- contrasts.fit(fit2, contrasts)
fit2C <- eBayes(fit2C)
topTable(fit2C, design2)
