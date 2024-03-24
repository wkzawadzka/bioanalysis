# Genomics

We share about 99.9% of our genetic makeup, yet small variations in our genomes create significant differences. These differences decide traits like height and influence disease susceptibility, such as cancer. It's the intricacies of our genomes that hold these answers. Let's learn more about genomics :)

## Week 1

### Mutation

Mutation is a change in your genome, often caused by DNA damage. During cell division, errors can occur in the replication process—about 1-3 errors per division. Occasionally, these errors can lead to harmful outcomes, like cancer, when they affect critical genes, disrupting their function.

### Central dogma

```
Information flows in signal direction:
from your genome (DNA) -> (transcription) -> to your DNA -> (translation) -> to proteins
```

- translation: read each triplet which is aminoacid (or 3 stop codon); 4^3 = 64 possible combinations

- proteins do most of the work in your body

- RNA: introns throw away, exons to proteins. So if you want to know what the proteins are that are being turned on in a cell, you need to know what the exons are.

### Cell biology

At most basic level, we divide cellular organisms into three domains:

- eukaryotes

- archaea

- bacteria

Archea and bacteria together called prokaryotes.
Eukaryotes have cell nuklei and prokaryotes have not.

One of the organelle's inside of the eukaryotic cell is called mitochondrion. Actually multiple in each cell. Has its own DNA. Powerhoue of the cell, as genes inside are resposible for a good bit of enery metabolism.
Fun fuct: The is a belief long long time ago mitochondrium was an independent prokaryote which was absorbed by the acnestor of all eukaryotic cells into the cell itself.

Cell division
(1) Mitosis - figure out how to divide copies

### Important molecules in molecular biology

- Adenine, Guanine -> purins, two ring structure

- Cutosine,Thymine -> one ring structure

A's bind to T's. G's bind to C's.

DNA in human genome is organized into 23 chromosome pairs. Each of these is a very, very long string. Longest in order of 250 mmillion nucleotides long!

DNA has a (1) direction (2) strandness.
One end of the DNA is 5' end and the other is 3' end. POsitive strand: 5' first, 3' second.

RNA "Template to make proteins" There are also "RNA genes"
-> no thymine, instead uracil.
A->A, G->G, C>C, T->Uracil; Single stranded

## Week 2

### PCR - Polymerase Chain Reaction

- very powerful way to copy DNA and make in fact, as many copies of DNA as you want to

**Primer**
-> a short sequence (like 15-20 bases) complementary to the sequence we wanna copy (so they stick-hybridize to it)

1. Heat up DNA (94C)
   It "melts" -> two strand fall apart, primer as well
2. Cool down gently (to 54C)
   Primer will stick to DNA first. Two strands would also come together BUT we don't wait that long.
3. (72C) Add copier molecule - DNA polymerase
   It looks for the parts of DNA where it is partly single stranded and partly double stranded and it grabs on to the sequence right there and starts to copy. Basically will start fill in missing places starting at the primer.

- Add lot's of A's, C's, G's and T's so polymerase has building blocks :)

### Next generation sequencing NGS

-> latest sequencing technology

2nd gen:
-> errors increase in later cycles (as we make million of copies in one of given "clusters" and we're sequencing that cluster one base at a time)
-> DNA polymerare adding a base isn't perfect
-> once in a while it will add a extra bse and given fragment will get ahead' or fail and get behind;

- exons sequencing: to know proteins

- RNA sequencing: to know what genes are being tuned up in a cell. We can't sequence RNA we need to turn in to DNA but we have very useful molecular mechanism: rerverse transcription.
  mrna

## Week 4

### "Scandal of Duke"

Case study: What went wrong in the Duke analysis?
Goal was: Try to decide which chemotherapies would apply best to which people on the besis of the genomic measurements.
(1) Lack of transparency
-> The data and code weren't reproducible
(2) Lack of cooperation
-> People in charge were reluctant to hand the data over to statisticians
== It takes a lot of time to discover problems in data analysis that are serious.
(3) Lack of expertise in statistics (they used silly prediciton rules)
(4) Major study design problems
-> batch effect (they run samples on different days)
-> different results in different days; reproducibility

### The Central Dogma of Statistics

**Knowing what the population is**

Say something about a population without measuring whole population by taking small sample of it.

- The best guess is not quite enough. How do we quantify possible variability?

### Data Sharing Plans

- 1. The raw dataset

- 2. A tidy data set (small, processed)

- 3. A code book describing each variable and its values in the tidy dataset

- 4. An expicit and exact recipe you used to go from 1 -> 2,3

** A guide for the lonely bioinformatician **
https://homolog.us/blogs/bioinfo/2013/04/23/todays-highlights-a-guide-for-the-lonely-bioinformatician/

### Plotting Your Data

```
Robert Gentleman, Genetech:
"Make big data as small as possible as quick as is possible" to enable sharing
```

- box plots: do plot data points

- plot replicates; be careful of scale - consider using transforms (...logs)

- Better: MA plots: rather than plotting one replicate on one axis and the other on another, you add the 2 duplicates od x-axis "Ave(log2R1, logR2)" and on y-axis you substract them. So on the left of the plot you have genes that are low expressed and on the right - that are very highly expressed. Looking at y-axis then, you look how far the replicates are from 0 and that signifies how different the two are from each other (so ideally if you actually measuring "replicates" you want it as close to 0 as possible)

Common phenomena in MA plots:
More differences between the replicates for the lowly expresed genes.

### Variability of genomic measurement

- phenotypic variability

- measurement error

- natural biological variation

### Statistical significance

When we know we observe difference between some groups large enough that we call it a "real" difference?

- t-statistic; average of group1
- average of group2 divided by a measure of variability. Denominator as scaling the difference between g1 and g2 by the units of variability.
  Small t-statisic: less likely that there's a difference
  Big t-statistic: more likely that there's a difference

- most used: p-value
  Idea: Suppose we calulated a t-statistic for comparing the difference between two groups. Suppose it equals to 2. Is it big value or small value?

-> permutations - randomly reassign labels and recaclulate the statistic; see the change. If you make the histogram of all the statistic you can get from random scrambles and you can see where the original statistic falls in that distibution.

P-value : sum up how many times the scrabled values were larger than the observed value (usually absolute value so you don't care about diffection of difference) and divide my # of permutations.

P-values closer to 0 are reported statistically significant.
Usuall cut-off 0.05.

**The p-value is not**

- Probability the null is true

- Probability that alternative is true

- A measure of statistical evidence

### Error rates

Suppose 550 out of 10,000 genes are significant at 0.05 level.

P-value < 0.05
Expect 0.05\*10,000 = 500 false positives
(P-values are uniform under the case where nothing's going on)
So we may think we found statistically significant results they might mostly be false positives.

Alternatively, suppose we used false discovery rate < 0.05.
Expect 0.05\*550 = 27.5 false positives

Finally, suppose we use family wise error rate < 0.05.
The probability of at least 1 false positive <= 0.05. That would mean almost all of them would be true positives.

Is P-value > 0.05 game over?

<img src="/" width="200" height="200">
