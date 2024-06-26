# course: Bioconductor for Genomic Data Science
# Taught by: Kasper Daniel Hansen, PhD
# Week 2 of Bioconductor for Genomic Data Science Johns Hopkins University

## DNA strings
library(Biostrings)
dna1 = DNAString("ACGT-G")
dna1
dna2=DNAStringSet(c("ACG", "ACGT", "ACGTT"))
dna2

# reversing DNA string vs DNA string set
rev(dna1) # seq: G-TGCA
rev(dna2) #    width seq
          #        5 ACGTT
          #[2]     4 ACGT
          #[3]     3 ACG


# translation
translate(dna1)

# freq
letterFrequency(dna2, letters = "GC")

library(BSgenome)
available.genomes()

BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer2")
library("BSgenome.Scerevisiae.UCSC.sacCer2")
Scerevisiae # genome of yeast
# Saccharomyces cerevisiae was the first eukaryotic genome
# to be completely sequenced.

seqnames(Scerevisiae)

seqlengths(Scerevisiae)

Scerevisiae$chrI

# letter frequency as percentage
letterFrequency(Scerevisiae$chrI, "GC", as.prob=TRUE)
# how to aggragate for all chromosomes? 
# apply function to all chromosomes by:
param = new("BSParams", X = Scerevisiae, FUN = letterFrequency)
# divide by sum of lengths of chromosomes
sum(unlist(bsapply(param, "GC"))) / sum(seqlengths(Scerevisiae))

# Now check a little bit with the GC content of the individual chromosomes.
unlist(bsapply(param, "GC", as.prob=TRUE))
      
### Biostrings - Matching- finding sequences & subsequences
dnaseq <- DNAString("ACGTACGT")

# single to single
matchPattern(dnaseq, Scerevisiae$chrI)
countPattern(dnaseq, Scerevisiae$chrI)

# 1 to many
vmatchPattern(dnaseq, Scerevisiae) # returned IRanges

# position weight matrix - search the genome for example for binding size for given transcription factor.
# matchPWM()

# Views
gr = vmatchPattern(dnaseq, Scerevisiae)
gr
vi = Views(Scerevisiae, gr)
vi # we get dna string :)

library(AnnotationHub)
ah = AnnotationHub()
qh = query(ah, c("sacCer2", "genes"))
qh
genes = qh[[1]]

# now lets get promotors
prom = promoters(genes)
prom = trim(prom)
promViews = Views(Scerevisiae, prom)
promViews


# "GC" content (nitrogenous bases)
gcProm = letterFrequency(promViews, "GC", as.prob=TRUE)
gcProm
plot(density(gcProm))


# GenomicRanges - Rle
# This is a form of compression of the very long vectors
# interesting for representing signal over the genome where the signal is only non-zero
# in a small part of the genome (ChIP seq...)

library(GenomicRanges)
rl = Rle(c(1,1,1,1,1,2,2,2,2,4,4,2,2))
rl
runLength(rl)

ir = IRanges(start = c(2,8), width = 4)
aggregate(rl, ir, FUN= mean)
vec = as.numeric(rl)
coverage(ir)


slice(rl, 2)
slice(rl, 3)


library(GenomicFeatures)

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
genes(txdb)

#subsetByOverlaps(genes(txdb), gr, ignore.strand=TRUE)
transcripts(txdb)
exons(txdb)
cds(txdb)


## rtracklayer - Data Import
# interfacing with a genome browser

library(rtracklayer)
ah = AnnotationHub()
table(ah$rdataclass)

ah.bw = subset(ah, rdataclass=="BigWigFile" & species == "Homo sapiens")
bw = ah.bw[[1]]
gr.ch22 = import(bw, which=GRanges("chr22", ranges = IRanges(1, 10^8)))


gr.ch22 = import(bw, which=GRanges("chr22", ranges = IRanges(1, 10^8)), as="Rle")


#What is the GC content of “chr22” in the “hg19” build of the human genome?
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)

chr22_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg19, "chr22")
nucleotide_freq <- letterFrequency(DNAString(chr22_sequence), letters = c("A", "G", "C", "T"))
(nucleotide_freq["G"] + nucleotide_freq["C"]) / sum(nucleotide_freq)


