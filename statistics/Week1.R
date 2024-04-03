#### Week1.R

tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19) # filled dots 

library(gplots)
library(devtools)
library(Biobase)
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)

con=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)

bm = bodymap.eset
pdata = pData(bm)
edata = exprs(bm)
fdata = fData(bm)

table(pdata$age)
# check if na
table(pdata$age, useNA="ifany")
# check if any missing values
sum(pdata$age==" ")
# check if na in samples
sum(is.na(edata)[1,])
# check in which row
rowSums(is.na(edata))
table(rowSums(is.na(edata)))
# check cols
table(colSums(is.na(edata)))

dim(fdata)
boxplot(log2(edata[,1]+1))
boxplot(log2(edata+1), col=2, range=0)

# disributions
plot(density(log2(edata[,1]+1)), col=2) # 1 sample
lines(density(log2(edata[,2]+1)), col=3)  # 2nd sample, diff color

#pltoDensities(log2(edata+1))

# qqplot
# check what sample has more values/
# less in given percentile

qqplot(log2(edata[,1]+1), log2(edata[,2]+1), col=2)
abline(c(0,1)) # 45 degree line


# difference between two samples
mm = log2(edata[,1]+1) - log2(edata[,3]+1)
aa = log2(edata[,1]+1) + log2(edata[,3]+1)
plot(aa, mm, col=2)
# if y=0 then no difference between samples
# each dot is one gene


### filter features (gene/probes) based on mean
edata = as.data.frame(edata)
filt_data <- filter(edata, rowMeans(edata)>1)
dim(filt_data)
boxplot(as.matrix(log2(filt_data+1)), col=2)



# ids for features
ids = as.character(fdata[,1])
ids[1:5]


chr = AnnotationDbi::select(org.Hs.eg.db, keys=ids, keytype = "ENSEMBL", columns="CHR")
head(chr)
dim(chr) # some got duplicated
chr = chr[!duplicated(chr[,1]),]
# check if match
all(chr[,1] == rownames(edata))
# take chromosome Y samples
edata <- as.data.frame(edata)
yChromData <- dplyr::filter(edata, chr$CHR=="Y")

boxplot(colSums(yChromData) ~ pdata$gender, col=2)
# overlay datapoints on top of ^ 
points(colSums(yChromData) ~ jitter(as.numeric(pdata$gender)),
       col=as.numeric(pdata$gender),
       pch=19)


# multivariate 
emat = as.matrix(edata)[rowMeans(edata)>10000, ]
heatmap(emat)

# blend the color from 2nd color from our tropical
# pallete, through white to 2nd color
colramp = colorRampPalette(c(3, "white", 2))(9)
heatmap(emat, col=colramp)

# without clustering
heatmap(emat, col=colramp, Rowv=NA, Colv=NA)

# with scale (how much differ?) - histogram color key
heatmap.2(emat, col=colramp, Rowv=NA, Colv=NA, dendrogram="none", scale="row", trace="none")

# ignore 0s
hist(log2(edata[,1]+1), breaks=100, col=2, xlim=c(1,15), ylim=c(0,400))

# for each gene, how many rows equal to 0?
hist(rowSums(edata==0),col=2)
# lots of genes that have all 0s
# few where all samples have non zero value

# we can remove lowly expressed genes 
low_genes = rowMeans(edata) < 5
table(low_genes)
# keep only those which are not low genes
filt_data = filter(as.data.frame(edata), !low_genes)

# by median
low_genes = rowMedians(edata) < 5
table(low_genes)
# keep only those which are not low genes
filt_data = filter(as.data.frame(edata), !low_genes)


#clustering
edata = data.frame(exprs(bm))
edata = edata[rowMeans(edata)>5000,]
edata = log2(edata + 1)
dist1 = dist(t(edata)) # distance between every set of samples

heatmap(as.matrix(dist1), col=colramp, Colv=NA, Rowv=NA)
hclust1 = hclust(dist1)
plot(hclust1, hang=-1)

# color dendograms
library(dendextend)
dend = as.dendrogram(hclust1)
dend = color_labels(hclust1, 4, 1:4)
plot(dend)

kmeans1 = kmeans(edata, centers=3)
matplot(t(kmeans1$centers), col=1:3, type="l", lwd=3)

# now plot heatmap in order according to cluster membership
newdata = as.matrix(edata)[order(kmeans1$cluster),]
heatmap(newdata,col=colramp, Colv=NA, Rowv=NA)
