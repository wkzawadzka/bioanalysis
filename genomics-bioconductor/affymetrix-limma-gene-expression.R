# course: Bioconductor for Genomic Data Science
# Taught by: Kasper Daniel Hansen, PhD
# limma affymatrix
library(limma)
BiocManager::install("leukemiasEset")
library(leukemiasEset)
data(leukemiasEset)

leukemiasEset # ann: genemapperhgu133plus2
table(leukemiasEset$LeukemiaType)
# ALL AML CLL CML NoL 
# 12  12  12  12  12 

data <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]
data$LeukemiaType # 5 levels
# let's correct levels
data$LeukemiaType <- factor(data$LeukemiaType)
data$LeukemiaType # 2 levels: ALL (leukemia) and NoL (no leukemia)


# design (how to set up-> see limma user's guide)
design <- model.matrix(~ data$LeukemiaType)
head(design) 
# intercept parameter is the average gene 
# expression for a given gene across ALL samples

# second column is the difference in gene expression
# between NOL and ALL
# if 0: gene is not diff expressed
# (no diff between NoL and ALL)

# --------- fit basic model for the design ---------
# h0: there is no difference between samples
# in the ALL group and the samples in the NOL group
# LM fits a linear model to all the genes separately

fit <- lmFit(data, design)
fit <- eBayes(fit) # esentially it shows that
# the variability of a gene is a mixture of a gene
# specific variability & global variability where var
# assumes that all genes have the same variability
# (not right biologically of course)
# but with only 12 samples, we do improve our variance
# estimation by having this mixture

topTable <- topTable(fit)

example <- topTable[1, ]

# let's compute by hand to check how it works:)
geneName <- rownames(example)[1]
# take one gene, do mean iside each level of the leukemia type
typeMean <- tapply(exprs(data)[geneName, ], data$LeukemiaType, mean)
typeMean # mean for ALL, mean for NoL
# log fold change
typeMean['NoL'] - typeMean['ALL']
# 4.089587  -> big


# --------- fit more specific model for the design ---------
design2 <- model.matrix(~ data$LeukemiaType - 1)
# new: -1
# why?
# it means we get new matrix where column names are
# slightly different 
head(design2)

colnames(design2) <- c("ALL", "NoL")
# cols:  data$LeukemiaTypeALL data$LeukemiaTypeNoL
# so expression level in ALL and expression level in NoL

# we need to make contrast
# contrast is hypothesis
# we want constast that is ALL - NoL
# this is opposite(!) than what we did before
# now: NoL is reference level.
contrast.matrix <- makeContrasts("ALL-NoL", levels=design2)
contrast.matrix

fit2 <- lmFit(data, design2)
fit2C <- contrasts.fit(fit2, contrast.matrix)
fit2C <- eBayes(fit2C)

topTable2 <- topTable(fit2C)
topTable2
# now its relative to our specific contrast 
# so we changed signe of our fold