library(BiocGenerics)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb) # The GenomeInfoDb package addresses a seemingly small, 
# but consistent problem: different online resources uses different
# naming conventions for chromosomes. In more technical jargon, thi
# package helps keeping your seqinfo and seqlevels harmonized. (Kasper Daniel Hansen, PhD)
# see: https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_seqinfo.html

############################## Genomic ranges ################################
### IRanges
ir <- IRanges(start = 1:3, width = 2)

# classic dataframe
df_ir_c <- data.frame(ir=ir)
df_ir_c

# DataFrame - more versatile
df_ir <- DataFrame(ir=ir, score = rnorm(3))
df_ir[1,1]


### GRanges
# ---
gr <- GRanges(seqnames= "chr1", strand = c('+', '-', '+'),
              ranges = IRanges(start = c(1,3,5), width = 3))
gr
values(gr) = DataFrame(score = rnorm(3))
values(gr)
gr$score
# ---
gr2 <- GRanges(seqnames= c("chr1", "chr2", "chr3"), strand = c("*"),
              ranges = IRanges(start = c(1,3,5), width = 3))
gr2
# ---
findOverlaps(gr, gr2)

# Very common use case: We have a GRange and we want to only select ranges
# from the GRange that overlaps some other elements:
subsetByOverlaps(gr, gr2)
# ---
### Seqinfo
# 'Even for analysis you may restrict you analysis to the six chromosomes
# or the autosomes. Or you may for convenience just want to look at chromosome
# 22, which is the smallest chromosome. So let's look at how we do that.
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
dropSeqlevels(gr, "chr2", pruning.mode="tidy")

############################### Annotation Hub #############################
# "AnnotationHub is an interface to a lot of different online resources. 
# The idea is that you create a hub, which is a local database of a lot of
# different online data that's out there. You take this local database,
# you query it, and you figure out which data do you want and then you go
# online and you retrieve them." (Kasper Daniel Hansen, PhD)
BiocManager::install("Bioconductor/AnnotationHub")
BiocManager::install("Bioconductor/BiocFileCache")

library(AnnotationHub)
ah = AnnotationHub()
ah[1]
unique(ah$dataprovider)


############################## Study ######################################
# - study of histone mark called H3K4 trimethylation, that is set to mark promoters

# 1. Take just human data
ah = subset(ah, species == "Homo sapiens")
# 2. Search for the hisotry modification data
# name of cell line: GM12878
# H3K4me3 is an epigenetic modification to the DNA packaging
# protein Histone H3 that indicates tri-methylation at the 4th
# lysine residue of the histone H3 protein
qhs = query(ah, c("H3K4me3", "Gm12878"))
qhs

gr1 = qhs[[2]] # broad peak 
gr2 = qhs[[4]] # narrow peak

summary(width(gr1))
summary(width(gr2))
table(width(gr2))

# Peaks are regions of the genome where there is an enrichment of ChIP-seq signal
peaks = gr2
# These peaks represent regions where a protein of interest (often a transcription factor or histone modification) is bound.
qhs[4]$genome # $genome: hg19

# Now, we're going to take these peaks, and we're going to ask, are these peaks enriched in promoters?

#  Promoters are regions of the genome that are upstream of genes and are
# involved in the initiation of transcription. When we say peaks are
# "enriched in promoters," it means we're investigating whether
# these regions of enrichment (peaks) are disproportionately located near
# gene promoters compared to what would be expected by chance.

query = query(ah, "RefSeq") # Reference Sequence Database
query$genome
genes = query[[1]]

# how many transcripts per gene
table(table(genes$name))


prom = promoters(genes)
prom
table(width(prom))

# Do we have some significant overlap?
# that is, is this particular histone modification, H3K4 trimethylation, enriched in promoters?
ov = findOverlaps(prom, peaks)
length(unique(queryHits(ov))) # 24826
length(unique(subjectHits(ov))) # 22073

# percentage of peaks that overlap promoter
length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE)) / length(peaks)

# 30% is pretty good but we need more analsis

length(subsetByOverlaps(prom, peaks, ignore.strand = TRUE)) / length(prom)
# 50%, that's pretty nice.


# Now, is this big or not? Let's look a little bit about the actual,
# let's look at how big are the peaks and how big are the promoters we
# have and how big is the intersection. 

sum(width(reduce(peaks, ignore.strand=TRUE))) ## 11 megabases
sum(width(reduce(prom, ignore.strand=TRUE))) ## 62 megabases
# And how big is the overlap?
sum(width(intersect(peaks, prom, ignore.strand=TRUE)))/10^6 # 3 megabases

inOut = matrix(0, ncol = 2, nrow = 2)
colnames(inOut) = c("in", "out")
rownames(inOut) = c("in", "out")
# both in peaks and promoters
inOut[1,1] <- sum(width(intersect(peaks, prom, ignore.strand=TRUE)))
# peaks on the rows, promoters on the columns
# a: in peaks, but not in promoters
inOut[1, 2] <- sum(width(setdiff(peaks, prom, ignore.strand = TRUE)))
# b: in promoters, but not in peaks
inOut[2, 1] <- sum(width(setdiff(prom, peaks, ignore.strand = TRUE)))
inOut
colSums(inOut)
rowSums(inOut)
# lets say human genome is 3m bases.(but is it okay? some of the bases have not even been sequenced)
inOut[2,2] <- 3*10^9 - sum(inOut)

# one diagonal by other
addsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
addsRatio # 18.26938
# Means the overlap is like 18 fold more enriched than we would expect

# Let's change our assumption. - We can do a sensitivity analysis
# So let's say that only half the human genome is mappable
inOut[2,2] <- 0
inOut[2,2] <- (3*10^9)/2 - sum(inOut) # much smaller genome
addsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
addsRatio # 8.915016
# Still much, much bigger than 1




############################## Test #################################
# Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.
ah = AnnotationHub()
ah = subset(ah, species == "Homo sapiens")
q = query(ah, "CpG")
q$genome
data_islands = q[[1]]

# How many islands exists on the autosomes?
# chr 1-22
autosome <- c(paste("chr", 1:22, sep=""))
autosomal_islands <- subset(data_islands, seqnames(data_islands) %in% autosome)
length(autosomal_islands)


# How many CpG Islands exists on chromosome 4?
autosomal_islands <- subset(data_islands, seqnames(data_islands) == "chr4")
length(autosomal_islands)

# Obtain the data for the H3K4me3 histone modification for the H1 cell
# line from Epigenomics Roadmap, using AnnotationHub. Subset these
# regions to only keep regions mapped to the autosomes (chromosomes 1 to 22).
# How many bases does these regions cover?
qhs = query(ah, c("H3K4me3", "E003"))
qhs
data_A <- qhs[[2]]
autosomal_A <- subset(data_A, seqnames(data_A) %in% autosome)
sum(width(autosomal_A))


# Obtain the data for the H3K27me3
qhs = query(ah, c("h3K27me3", "E003"))
qhs
data <- qhs[[2]]
autosomal_B <- subset(data, seqnames(data) %in% autosome)
#  In the return data, each region has an associated "signalValue". 
# What is the mean signalValue across all regions on the standard chromosomes?
mean(autosomal_B$signalValue)


# Bivalent regions are bound by both H3K4me3 and H3K27me3.
# Question: Using the regions we have obtained above,
# how many bases on the standard chromosomes are bivalently marked?
sum(width(intersect(autosomal_A, autosomal_B)))

  
# Question 6
# We will examine the extent to which bivalent regions overlap CpG Islands.
# Question: how big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap one or more CpG Islands?
bivalent = intersect(autosomal_A, autosomal_B)
length(subsetByOverlaps(bivalent, autosomal_islands, ignore.strand = TRUE)) / length(bivalent)


# Question 7
# Question: How big a fraction (expressed as a number between 0 and 1)
# of the bases which are part of CpG Islands, are also bivalent marked.
biv_islands = intersect(bivalent, autosomal_islands)
sum(width(biv_islands)) / sum(width(autosomal_islands))

# 8 Question: 
# How many bases are bivalently marked within 10kb of CpG Islands?
# Tip: consider using the "resize()"" function.
autosomal_inslands_10kb <- resize(autosomal_islands, width = width(autosomal_islands) + 20000, fix="center")
sum(width(intersect(autosomal_inslands_10kb, bivalent)))

  
  