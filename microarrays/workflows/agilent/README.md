# Dual probe Agilent workflow

### About Agilent microarrays

Agilent combines two-sample hybridization with the use
of `long (60-mer) oligonucleotides`. [1] Which is different from Affymetrix short-length (about 25 oligonucleotides) probes.

### Dual probe

In the two-color cDNA platform, mRNA from two samples is reverse-transcribed to cDNA, labeled with fluorescent dye, `Cy3 (green)` - usualy for control group and `Cy5 (red)`, and simultaneously
hybridized to an array containing spots of DNA
sequences. `Ratios of fluorescence intensities` provide a relative measure of expression at each spot on the array. [1]

### Datasets

For now, 3 Agilent datasets are used in the analysis: GSE11223, GSE179285, GSE96665. GSE111761 was not used as it only have samples with CD, not UC.

## Pre-processing

Goal: retaining biological variation while minimizing artificial systematic (labeling bias, irregular feature morphologies, mismatched sample concentrations, cross-hybridization [1, 4]) distortions.

### Read raw intensities - foreground estimators

The intensity of each spot is summarized by median pixel intensity of the spot. Experiments [1] shown the choice between mean and median had little effect on the final analysis.

`RG <- read.maimages(
  targets,
  source = 'agilent'
)`

The use of `source="agilent"` defaults to `"agilent.median"`

### Normalization

**(A) Within-array normalization**

Considering each sample (array) separately.

**(B) Between-array normalization**

### Background correction

Background correction attempts to `adjust the data for the ambient intensity surrounding each feature`. The “normexp” method models the observed pixel intensities as the sum of 2 random variables, one normally distributed and the other exponentially distributed, representing background noise and signal, respectively. [2] Using a saddle-point approximation, Ritchie and others [3] found `normexp to be the best background correction method for 2-color microarray data`.

### Adjusting for False Discovery Rate (FDR)

To avoid getting FP we need to correct p-values, as the number of featues for which testing whether the gene is significantly differentially expressive is 30,000+, which would mean (30000\*0.0.5) 1500 are DEGs by chance. Method used - Benjamini-Hochberg.

### Results

1. **Downregulated DEGs**

| Gene     | LogFC  | adj.P.val | Description                                                                                                                                                                               | Additional info                                                                                              |
| -------- | ------ | --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| LGALS2   | -0.757 | 4.621e-03 | oxidative stress-responsive gene with an inhibitory function on colon tumor growth                                                                                                        | [A](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7790754/) \| [B](https://pubmed.ncbi.nlm.nih.gov/25619132/) |
| CLDN8    | -0.701 | 0.019     | partially determine the integrity and paracellular permeability of the intestinal epithelium; <br> downregulation of claudin 8 has been observed in intestinal inflammatory disorders [B] | [A](https://pubmed.ncbi.nlm.nih.gov/23205909/) \| [B](https://pubmed.ncbi.nlm.nih.gov/28493289/)             |
| BC032913 | -0.901 | 5.441e-05 | antisense lncRNA; is Downregulated in Gastric Cancer [A]                                                                                                                                  | [A](https://pubmed.ncbi.nlm.nih.gov/32914372/)                                                               |

### Citations

[1] Zahurak, M., Parmigiani, G., Yu, W. et al. Pre-processing Agilent microarray data. BMC Bioinformatics 8, 142 (2007). https://doi.org/10.1186/1471-2105-8-142

[2] Silver JD, Ritchie ME, Smyth GK. Microarray background correction: maximum likelihood estimation for the normal-exponential convolution. Biostatistics. 2009 Apr;10(2):352-63. doi: 10.1093/biostatistics/kxn042. Epub 2008 Dec 8. PMID: 19068485; PMCID: PMC2648902.

[3] Ritchie ME, Silver J, Oshlack A, Holmes M, Diyagama D, Holloway A, Smyth GK. A comparison of background correction methods for two-colour microarrays. Bioinformatics. 2007 Oct 15;23(20):2700-7. doi: 10.1093/bioinformatics/btm412. Epub 2007 Aug 25. PMID: 17720982.

[4] Delenstarr G, Cattel H, Chen C, Dorsel A, Kincaid R, Nguyen K, Sampas N, Schidel S, Shannon K, Tu A, Wolber P: Estimation ofthe
confidence limits of oligo nucleotide array-based measurements of differential expression. SPIE Proceedings: Microarrays:
Optical TEchnologies and Informatics 4266 2001:120-131.
