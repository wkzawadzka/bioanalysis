library(devtools)
library(Biobase)
library(ggplot2)
library(broom)
library(snpStats)
library(MASS)
library(DESeq2)

tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
pch(19)

data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 = snps.10[,use]
## xx transpose
xxmat <- xxt(sub.10, correct.for.missing = FALSE)
# eigen vectors
evv <- eigen(xxmat, symmetric = TRUE)
pcs <- evv$vectors[,1:5]
dim(pcs)


snpdata <- sub.10@.Data
status <- subject.support$cc
snp1 <- as.numeric(snpdata[,1])
table(snp1)
snp1[snp1=0] = NA

# fit model
glm1 <- glm(status ~ snp1, family="binomial")
glm1

# dominant model
snp1_dom = (snp1==1)

glm1 <- glm(status ~ snp1_dom, family="binomial")
glm1

# adjusting for principal components
glm2 <- glm(status ~ snp1 + pcs[,1:5], family="binomial")
tidy(glm2)


glm_all <- snp.rhs.tests(status ~ 1, snp.data=sub.10)
slotNames(glm_all)
qq.chisq(chi.squared(glm_all), df=1) # qq 
# we see the line is not great so let's adjust for pcs
glm_all <- snp.rhs.tests(status ~ pcs, snp.data=sub.10)
slotNames(glm_all)
qq.chisq(chi.squared(glm_all), df=1) # qq nice :)


# quantifying relationships
library(genefilter)
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
save(bottomly.eset, file="bottomly.Rdata")

load(file="bottomly.Rdata")
ls()
bat <- bottomly.eset
pdata <- pData(bat)
edata <- as.matrix(exprs(bat))
fdata <- fData(bat)
edata <- log2(as.matrix(edata)+1)
edata <- edata[rowMeans(edata)>10,]

# t statistic
# all genes at once with genefilter package
tstats_obj <- rowttests(edata, pdata$strain)
names(tstats_obj)
hist(tstats_obj$statistic, col=2)

fstats_bj <- rowttests(edata, as.factor(pdata$lane.number))


# moderated t statistic
mod <- model.matrix(~ pdata$strain)
fit <- lmFit(edata, mod)
fit <- eBayes(fit)
head(fit$t)


# adjusted
mod_adj <- model.matrix(~ pdata$strain  + as.factor(pdata$lane.number))
fit2 <- lmFit(edata, mod_adj)
fit2 <- eBayes(fit2)
head(fit2$t)

# permutation
hist(tstats_obj$statistic, col=1, xlim=c(-5,2)) # previous
set.seed(222)
strain0 <- sample(pdata$strain)
tstats_obj0 <- rowttests(edata, strain0)
hist(tstats_obj0$statistic, col=1, xlim=c(-5,2))
# new statistics are mostly positive (often because there is
# a covariate that is not modelled)
quantile(tstats_obj0$statistic)
quantile(tstats_obj$statistic) 



library(limma)
BiocManager::install('edge')
library(edge)
library(qvalue)


fstats <- rowFtests(edata, as.factor(pdata$strain))
hist(fstats$p.value, col=1) # we see it looks wrong. we need better model

# another method
edge.study <- build_study(edata, grp = pdata$strain, adj.var = as.factor(pdata$lane.number))
diff_expr <- lrt(edge.study)
qval <- qvalueObj(diff_expr)
hist(qval$pvalues, col=4)

# modearte statistic
mod <- model.matrix( ~ pdata$strain + pdata$lane.number)
fit <- lmFit(edata, mod)
fit <- eBayes((fit))
pvals <- topTable(fit, number=Inf)$P.Value
hist(pvals, col=3)


