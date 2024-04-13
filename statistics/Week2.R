### Week 2 Statistics for Genomic Data Science
# Johns Hopkins University Taught by: Jeff Leek, PhD

library(devtools)
library(Biobase)

par(pch=19) # filled dots 
tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')

palette(tropical)


con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)

mp = montpick.eset
pdata = pData(mp)
edata = as.data.frame(exprs(mp))
fdata = fData(mp)

# reduce size
edata <- edata[rowMeans(edata)>100, ]
dim(edata) # 3072  129
edata <- log2(edata+1)

### SVD
# We need to center the data as we will be doing SVD
# If you dont remove it then the first singular value vector will always be mean level
# since that will always explain most deviation in genomic data
edata_centered <- edata - rowMeans(edata) # genes
svd1 = svd(edata_centered)
names(svd1)
dim(svd1$u) #  3072  129 # variation across samples
length(svd1$d) # 129
dim(svd1$v) # 129 129 # variation across genes


plot(svd1$d, ylab="Singular values", col=2)
# plot the variance explained (percentge)
plot(svd1$d^2/sum(svd1$d^2), ylab="Percent Variance Explained", col=2)
# We can see ^ that the first singular value explains more than 50% of the variance
# (highly explanatory variable) Let's see:
par(mfrow=c(1,2))
plot(svd1$v[,1], col=1, ylab="1st PC")
plot(svd1$v[,2], col=3, ylab="2nd PC")

# plotting first signular vector against second singular vector
par(mfrow=c(1,1))
# col -> color by what study they come from
plot(svd1$v[,1], svd1$v[,2], ylab="2nd PC", xlab="1st PC", col=as.numeric(pdata$study))
# we can also do the boxplot
boxplot(svd1$v[,1] ~ pdata$study, border=c(1,2))
# overlay points over boxplot
points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)), col=as.numeric((pdata$study)))

### Principal Component Analysis
pc1 <- prcomp(edata) # ! columns
plot(pc1$rotation[,1], svd1$v[,1]) # not the same thing as not scaled in the same way
# lets substact col means now
edata_centered2 = t(t(edata)-colMeans(edata)) # centred by column

svd2 = svd(edata_centered2)
plot(pc1$rotation[,1], svd2$v[,1]) # not they are identical


### Outlier analysis
edata_outlier <- edata_centered
edata_outlier[6,] = edata_outlier[6, ]*1000 # make an outlier out of 6th gene
svd3 <- svd(edata_outlier)

plot(svd1$v[,1], svd3$v[,1], xlab="Without outlier", ylab="With outlier") # they dont match each other anymore
# no match in terms of their signular value decomposition
# let's see how first component reflects the outlier
plot(svd3$v[,1], edata_outlier[6,], col=3) # highly correlated!
# What is happening is that the decomposition: if one gene expression is way higher
# than the others are, it is gonna drive most the variation in the dataset
# CONCLUSION: Take care of scaling & centering before SVD.

### Quantile normalization
library(preprocessCore)

edata = as.data.frame(exprs(mp))
edata = log2(edata+1)
edata = edata[rowMeans(edata)>3, ]
dim(edata)

colramp = colorRampPalette(c(3, "white", 2))(20)
plot(density(edata[,1]), col=colramp[1], lwd=3, ylim=c(0,.30))
for (i in 2:20){lines(density(edata[,i]), lwd=2, col=colramp[i])}
# if differences likely to technology not biology
norm_data <- normalize.quantiles(as.matrix(edata))
dim(norm_data)
for (i in 1:20){lines(density(norm_data[,i]), lwd=2, col=colramp[i])}

plot(norm_data[1,], col=as.numeric(pdata$study))
svd1 = svd(norm_data - rowMeans(norm_data))
plot(svd1$v[,1], svd1$v[,2], col=as.numeric(pdata$study)) # variability not removed


### Linear models
library(broom)
# Illumina BodyMap data
con=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)

bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()

# (1) numerical data - age
edata = as.matrix(edata)
lm1 = lm(edata[1,] ~ pdata$age) # first gene
tidy(lm1)
# plot age vs first gene
plot(pdata$age, edata[1,], col=1)
abline(lm1, col=2, lwd=3)

# (2) categorical data - sex
table(pdata$gender)
boxplot(edata[1,] ~ pdata$gender)
points(edata[1,] ~ jitter(as.numeric(pdata$gender)), col=as.numeric(pdata$gender))

dummy_m = pdata$gender == 'M'
dummy_m * 1

lm2 <- lm(edata[1,] ~ pdata$gender)
tidy(lm2)

lm3 = lm(edata[1,] ~ pdata$age + pdata$gender)
tidy(lm3)
summary(lm3)

# age interacting with gender
# tricky in interpreting
lm4 = lm(edata[1,] ~ pdata$age*pdata$gender)
tidy(lm4)
summary(lm4)

lm5 = lm(edata[6,] ~ pdata$age)
plot(pdata$age, edata[6,], col=4)
abline(lm5, col=5, lwd=3)
# outlier did not affect the line

index <- 1:19
lm6 <- lm(edata[6,] ~ index)
plot(index, edata[6,]) # now outlier at the end
abline(lm6, col=1, lwd=3)
# line affected
# lets see without this data point
lm7 <- lm(edata[6,-19] ~ index[-19])
abline(lm7, col=3, lwd=3)
legend(5, 1000, c("With outlier", "without outlier"), col=c(1,3), lwd=3)

# let's make histogram of the residuals
par(mfrow=c(1,2))
hist(lm6$residuals, col=2) # huge outlier residuals
hist(lm7$residuals, col=3)

gene1 = log2(edata[1,]+1)
lm8 <- lm(gene1 ~ index)
hist(lm7$residuals) # look less bad


model.matrix( ~ index)
dim(model.matrix( ~ index))

colramp <- colorRampPalette(1:4)(17)
lm9 <- lm(edata[2,] ~ pdata$age)
par(mfrow=c(1,1))
plot(lm9$residuals, col=colramp[as.numeric(pdata$tissue.type)])
legend(10,2000, pdata$tissue.type, col=colramp[as.numeric(pdata$tissue.type)], lwd=5)

# Many Regressions at Once
# https://jtleek.com/genstats_site/lecture_notes/02_10_Many_Regressions.pdf
# Adjustment variables: where data comes from,
# gender so on

# Primary (“Biological”) Variable
# Adjustment Variables

# Obvious:
# Location (Agadir, Desert, Village)
# Sex (female, male)
# Batch (data generated in two batches)
# Subtle:
# Intensity dependent effects
# Dye effects
# Probe composition effect


# Many Regressions
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()


edata <- log2(as.matrix(edata)+1)
edata <- edata[rowMeans(edata)>10,] # remove lowly expressed genes
dim(edata)

mod <- model.matrix( ~ pdata$strain)
fit <- lm.fit(mod, t(edata)) # fit for every single gene
names(fit)
fit$coefficients[,1]

par(mfrow=c(1,2))
hist(fit$coefficients[1,], breaks=100, col=2, xlab='Intercept')
hist(fit$coefficients[2,], breaks=100, col=2, xlab='Strain')

# adjusting for lane
mod2 = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_adj <- lm.fit(mod2, t(edata))

library(limma)
fit_limma = limma::lmFit(edata, mod2)
names(fit_limma)



### Batch effects  and Confounders
# eBayes for adjusting for batch effects
BiocManager::install('snpStats')
library(sva)
library(bladderbatch)
library(snpStats)

data(bladderdata)
ls()
pheno = pData(bladderEset)
edata = as.matrix(exprs(bladderEset))
head(pheno) # cel files

# Method 1: just adjust in the model
mod <- model.matrix(~ as.factor(pheno$cancer) + as.factor(pheno$batch))
fit <- lm.fit(mod, t(edata))
par(mfrow=c(1,1))
hist(fit$coefficients[2,], col=1, breaks=100)
# be careful of the correlation of cancer & batch
table(pheno$cancer, pheno$batch)
# ideally: same sizes - but its not the case.

# Method 2: ComBat
# a model matrix for combat (just intercept term)
modcombat <- model.matrix( ~1, data=pheno)
# b model we are gonna test
modcancer <- model.matrix(~pheno$cancer, data=pheno)

combat_edata <- ComBat(dat=edata, batch=pheno$batch, mod=modcombat, par.prior=T, prior.plots = TRUE)
combat_fit <- lm.fit(modcancer, t(combat_edata))
hist(combat_fit$coefficients[2,], col=2, breaks=100)

plot(fit$coefficients[2,], combat_fit$coefficients[2,], col=2, xlim=c(-5,5), ylim=c(-5,5))
abline(c(0,1), col=3, lwd=3)
# combat values smaller so x axis smaller than y

### SVA::  Case: we dont know batch effects
mod <- model.matrix(~ pheno$cancer, data=pheno) # model including the variable we care about
mod0 <- model.matrix(~1, data=pheno)

sval <- sva(edata, mod, mod0, n.sv=2)
names(sval)
summary(lm(sval$sv ~ pheno$batch))
# second one super correlated with batch
# pheno$batch  -3.994 0.000194 ***

modsv <- cbind(mod, sval$sv)
head(modsv)
fitsv <- lm.fit(modsv, t(edata))
par(mfrow=c(1,2))
# compare sva to combat
plot(fitsv$coefficients[2,], combat_fit$coefficients[2,], col=4, xlim=c(-5,5), ylim=c(-5,5))
abline(c(0,1), lwd=3, col=3)
# compare sva to adjusted linear model
plot(fitsv$coefficients[2,], fit$coefficients[2,], col=4, xlim=c(-5,5), ylim=c(-5,5))
abline(c(0,1), lwd=3, col=3)
abline(c(0,1), lwd=3, col=3)

### Genetic dana, different populations
data(for.exercise)
ls()
controls <- rownames(subject.support)[subject.support$cc==0]
use <- seq(1, ncol(snps.10), 10) #every 10th values
ctl.10 <- snps.10[controls, use]

# calulate principal components
xxmat <- xxt(ctl.10, correct.for.missing = FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5] # first 5
dim(pcs)

# look what population they come from
pop <- subject.support[controls, subject.support$stratum] # population stratification variable
head(pop)
# plot first PC to second PC
par(mfrow=c(1,1))
plot(pcs[,1], pcs[,2], col=as.numeric(pop$stratum), xlab="PC1", ylab="PC2")
legend(-0.02,0.15, legend=levels(pop$stratum), pch=19, col=1:2)



### QUIZ
# 1. Question 1 Load the Montgomery and Pickrell eSet:
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
# What percentage of variation is explained
# by the 1st principal component in the data set if you:
#   Do no transformations?
#   
#   log2(data + 1) transform?
#   
#   log2(data + 1) transform and subtract row means?

a_data = edata
b_data = log2(edata+1)
c_data = log2(edata+1) - rowMeans(log2(edata +1))

pc_a <- prcomp(a_data)
pc_b <- prcomp(b_data)
pc_c <- prcomp(c_data)

pc_a$sdev^2/sum(pc_a$sdev^2)
pc_b$sdev^2/sum(pc_b$sdev^2)
pc_c$sdev^2/sum(pc_c$sdev^2)

# 2. Question 2
# Perform the log2(data + 1) transform and subtract row means from the samples. Set the seed to 
# 333
# 333 and use k-means to cluster the samples into two clusters. Use 
# svd
# svd to calculate the singular vectors. What is the correlation between the first singular vector and the sample clustering indicator?

logged_data <- log2(edata+1)
centered_data <- logged_data - rowMeans(logged_data)
set.seed(333)
svd <- svd(centered_data)
km <- kmeans(t(centered_data), centers = 2)
cor.test(svd$v[,1], km$cluster)

# 3. Question 3 Load the Bodymap data with the following command
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

.
Question 3
Load the Bodymap data with the following command

3456
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
# Fit a linear model relating the first gene’s counts to the number of technical replicates, treating the number of replicates as a factor. Plot the data for this gene versus the covariate. Can you think of why this model might not fit well?

lm <- lm(edata[1,] ~ as.factor(pdata_bm$num.tech.reps))
plot(pdata_bm$num.tech.reps, edata[1,], col=3)
#abline(lm$coefficients, col=2, lwd=3)
summary(lm)
tidy(lm)


# 4. Question 4
Load the Bodymap data with the following command

456
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
# Fit a linear model relating he first gene’s counts to the age of the person and the sex of the samples. What is the value and interpretation of the coefficient for age?

lm2 <- lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
tidy(lm2)
lm2$coefficients


# 5 Qs
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)

edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata <- log2(edata+1)
mod = model.matrix(~ pdata$population, data=pdata)

fit <- lm.fit(mod, t(edata))
dim(fit$residuals)
dim(fit$effects)
dim(fit$coefficients)

# 6. Q
?lm.fit
# n vector of orthogonal single-df effects. The first rank of them correspond to non-aliased coefficients, and are named accordingly.

# 7.
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

isAgeNull <- is.na(pdata_bm$age)
sum(isAgeNull)
mod = model.matrix(~ pdata_bm$age[!isAgeNull])
lm <- limma::lmFit(edata[,!isAgeNull], mod)
lm$coefficients[1000,]
plot(edata[1000,])
abline(lm)
summary(lm)

# 8.
# # Fit many regression models to the expression data where 
# age
# age is the outcome variable and 
# tissue.type
# tissue.type is an adjustment variable using the 
# lmFit
# lmFit function from the 
# limma
# limma package (hint: you may have to subset the expression data to the samples without missing values of age to get the model to fit). What is wrong with this model?

mod <- model.matrix(~ pdata_bm$age[!isAgeNull] + as.factor(pdata_bm$tissue.type[!isAgeNull]))
fits <- limma::lmFit(edata[,!isAgeNull], mod)


# 10
set.seed(33353)
edata = exprs(bm)
edata <- as.matrix(edata)
edata <- log2(edata +1)
edata <- edata[rowMeans(edata)>1, drop=FALSE]
dim(pdata_bm)
mod = model.matrix(~ pdata_bm$age[!isAgeNull], data=pdata_bm)
mod0 <- model.matrix(~1, data=pdata_bm)

sva1 <- sva(edata, mod, mod0, n.sv=2)
cor(sva1$sv, pdata_bm$age)