## --------------------------------------------------------------
accesion_number <- "GSE11223"
setwd("C:/Users/weraz/Pictures/Sem6/R-bioconductor/tests/data/agilent")

library(GEOquery)
library(tidyverse)
library(data.table)
library(ggplot2)
library(limma)
library(biomaRt)
library(arrayQualityMetrics)
library(testit)
library("gplots")
library("AgiMicroRna")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "curl")
options

tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

# isLog2Transformed <- function(data){
#   qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#   shouldBeLogged <- (qx[5] > 100) ||
#     (qx[6]-qx[1] > 50 && qx[2] > 0)
#   return(!shouldBeLogged)
# }

## download ----------------------------------------------------
#getGEOSuppFiles(accesion_number)
#tar_file <- paste0(accesion_number, "/", accesion_number, "_RAW.tar")
#untar(tar_file, exdir = paste0(accesion_number, "/raw/"))

## get file names  ----------------------------------------------
file_paths <- list.files(file.path(accesion_number, "raw"), pattern = "GSM.*\\.gz", full.names = TRUE)
cat("There are", length(file_paths), "individual biological samples in the", accesion_number, "database.\n")

# create targets ------------------------------------------------
df_pheno <- fread(paste0("./", accesion_number, "/", accesion_number, "_newPheno.txt.gz"))
assert(tools::file_path_sans_ext(tools::file_path_sans_ext(basename(file_paths))) == df_pheno$GSM)
df_target <- data.frame(
  FileName = file_paths,
  Treatment = ifelse(df_pheno$Category == 0, "Control", "UC"), 
  GErep = ifelse(df_pheno$Category == 0, 1, 2)
)
# read ---------------------------------------------------------
RG <- read.maimages(
  df_target,
  source="agilent"
) 
# exploratory analysis -----------------------------------------
table(rowSums(is.na(RG$G)))
table(colSums(is.na(RG$G)))
gene_freq = as.data.frame(table(RG$genes$GeneName, useNA="ifany"))
gene_freq = gene_freq[order(-gene_freq$Freq),]
gene_freq[1:20,]
ggplot(gene_freq[1:20,], aes(x = Freq, y = reorder(Var1, Freq))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = Freq), hjust = -0.1, size = 3) +
  labs(x = "Frequency", y = "Probe/Gene name") +
  theme_minimal() # most control 

# "3 tables"
# check if # of rows in pdata is same and cols in edata
# and rows of edata match rows of fdata
assert(dim(RG$G)[2]==dim(df_pheno)[1])
assert(dim(RG$G)[1]==dim(RG$genes)[1])
# annotate prior to DEG analysis -------------------------------
# prior as we want to save our processing to file :)
annotations <- fread(file.path(accesion_number, "raw", "GPL1708_old_annotations.txt.gz"))
annotations <- annotations[which(annotations$SPOT_ID %in% RG$genes$ProbeName),]
annotations <- annotations[match(RG$genes$ProbeName, annotations$SPOT_ID),]
table(RG$genes$ProbeName == annotations$SPOT_ID) # check allignment
RG$genes$Entrez <- annotations$GENE
#RG$genes$Ensembl <- annotations$ENSEMBL_ID
RG$genes$GeneName <- annotations$GENE_SYMBOL
RG$genes$Name <- annotations$NAME
RG$genes$ID <- annotations$ID
# preprocess ---------------------------------------------------
RG.b <- backgroundCorrect(RG,  method="normexp", offset = 20)
MA <- normalizeWithinArrays(RG.b, method="loess")
MA.n <- normalizeBetweenArrays(MA, method="Aquantile")
#MA.avg <- avereps(MA.n, ID=MA$genes$ProbeName)
MA.n <- MA.n[!(is.na(MA.n$genes$Name)), ]
# genes with multiple "copies" across the genome 
# todo: averaging or take one with best FC?
#MA.avg <- avereps(MA.n, ID=MA.n$genes$Name) # Probe name from annotations
#dim(MA.n) # 41689   202
#dim(MA.avg) # 33492   202 

# account for multiple probes
MA.avg = avereps(MA.n, ID=MA.n$genes$ProbeName)
dim(MA.n) # 41689   202
dim(MA.avg) # 41001   202

# multivariate sample comparison across genes ------------------
colramp = colorRampPalette(c(3, "white", 2))(9)
procedureDates <- as.Date(df_pheno$Procedure.Date, format = "%m/%d/%y")
MOrdered <- MA.n$M[,order(procedureDates)]
procedureDates[29:30] # same: "2004-12-02" "2004-12-02"
heatmap(MOrdered[1:20000,28:31], col=colramp, Rowv=NA, Colv=NA)
MA.n$targets$FileName[order(procedureDates)][29:30] # GSM282877 GSM282891
MA.n$targets$Treatment[order(procedureDates)][29:30] # "Control" "Control"
# GSM282891	11376 Normal Uninflamed ascending colon
# GSM282877	11082 Normal Uninflamed descending colon
# possible resason:
# the most expressed one is from right colon, all rest from left.

kmeans1 = kmeans(MA.n$M, centers=6)
matplot(t(kmeans1$centers), col=1:5, type="l", lwd=3)
newdata = as.matrix(MA.n$M)[order(kmeans1$cluster),]
#heatmap(newdata,col=colramp, Colv=NA, Rowv=NA) #too long


# preprocessing comparison of intensity densities --------------
par(mfrow = c(2, 2), mar = c(3, 2, 2, 2))
title_size <- 0.8
plotDensities(RG, main = "", legend = FALSE)
title(main = "RG Densities - raw", cex.main = title_size)
plotDensities(RG.b, main = "", legend = FALSE)
title(main = "Background corrected", cex.main = title_size)
plotDensities(MA, main = "", legend = FALSE)
title(main = "Within array normalized", cex.main = title_size)
plotDensities(MA.n, main = "", legend = FALSE)
title(main = "Within- and between-\narray normalized", cex.main = title_size)
par(mfrow = c(1, 1))
# save matrix --------------------------------------------------
# not sure about this 
# todo: what is geo series matrix rows selection?
Control <- MA.avg$genes$ControlType==1L
MA.filtered <- MA.avg[!Control,  ]
dim(MA.filtered)
dim(MA.avg)
#df_pheno$GSM[1:5]
write.table(MA.filtered$M, file = paste0(accesion_number, "/matrix.txt"), row.names = MA.filtered$genes$ID, col.names = df_pheno$GSM)
# gene expression ----------------------------------------------
design <- model.matrix(~ df_target$Treatment - 1)
colnames(design) <- c("Control", "UC")
cont.matrix <- makeContrasts('UC-Control', levels = design)
fit <- lmFit(
  MA.filtered,
  coef = "UC-Control",
  design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
tt <- topTable(fit, adjust = 'BH', number = Inf)
# A value of p < 0.01 and a fold change of greater or less than 1.5 were considered statistically significant.
# & Subset significant genes based on FDR threshold (paper: 5%)
tt_subset <- subset(tt, adj.P.Val < 0.05 & P.Value < 0.01 & abs(logFC)>log2(1.5))
tt_subset <- tt_subset[, c("Entrez","GeneName", "Description", "logFC", "adj.P.Val", "P.Value")]

# just with entrez
tt_entrez <- tt_subset[!is.na(tt_subset$Entrez),]
#row.names(tt_entrez) <- tt_entrez$Entrez
tt_entrez <- tt_entrez[,c("Entrez", "GeneName", "logFC", "adj.P.Val", "P.Value", "Description")]
# plots ------------------------------------------------------
hist(tt$logFC, breaks = 50, main = "Distribution of Fold Changes", xlab = "Log2 Fold Change")
hist(x=fit$p.value, freq=TRUE, main="Histogram of p values")
# decide test ------------------------------------------------
# alternatively:
# results <- decideTests(fit, method="global", lfc=log2(1.5), p.value = 0.05)
# summary(results)
# vennDiagram(results, circle.col=palette(tropical))
# probeNames <- rownames(results[which(results[,1]!=0),])
# probeNames
# tt_genes <- subset(tt, SystematicName %in% probeNames | ProbeName %in% probeNames)[, c("Description", "GeneName", "logFC", "adj.P.Val", "P.Value")]
# ------------------------------------------------------------
write.fit(fit, adjust="BH", method="global", file=paste0(accesion_number, "/globalDEresults.txt"))
# ------------------------------------------------------------
par(mfrow = c(3,5))
par(mgp = c(2, 0.4, 0))
par(mar = c(4, 4, 1, 1) + 0.1)
for (column in 150:length(colnames(MA))) {
  group <- df_target$Treatment[column]
  name <- df_target$Label[column]
  #group <- targets$Treatment[grep(paste0("^", filename), targets$FileName)]
  #name <- gsub("\\.txt$", "", basename(filename))
  
  limma::plotMD(MA, column = column, 
                xlab = "Avg log-expression\n(means)", 
                ylab = "Expression log-ratio\n(this sample vs others)", 
                main = paste0(name, " (", group, ")"),
                cex.main = 0.7,
                cex.lab = 0.8, 
                cex.axis = 0.8  # tick
  )   
}
par(mfrow = c(1, 1))
