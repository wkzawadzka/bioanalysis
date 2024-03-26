# Bioconductor for Genomic Data Science

## Course by Kasper Daniel Hansen, PhD

## Week 1

### GRanges

**GRanges** is a data structure for storing genomic intervals.
IRanges is a vector that contains integer intervals. GRanges similar but having to do with chromosome and strength.

Rest of the notes in `Week1.R` - GRanges, AnnotationHub.

- study of histone mark called H3K4 trimethylation, that is set to mark promoters

## Week 2

`Week2.R`

- Biostrings

- BSgenome

- GenomicRanges - Rle (form of compression of the very long vectors); interesting for representing signal over the genome where the signal is only non-zero in a small part of the genome;

- GenomicRanges - Lists

- GenomicFeatures

- rtracklayer - Data Import; interfacing with a genome browser;

## Week 3

Annotation: linking our experimental data to various databases or repositories.

- Entrez
  ENTREZ is a metadatabase hosted by NCPI which is kind of a union of a set, or a superset of a lot of different databases that are hosted by the federal government.

Two main approaches to annotation in Bioconductor: package ex.hgu95av2.db or query online resources like UCSC or ENSEMBL (latest version through online resources!).

- Expression matrix
  Convention: samples on columns and features (genes or probes) on the rows
