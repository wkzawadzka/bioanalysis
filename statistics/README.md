# Statistics for Genomic Data Science

# Johns Hopkins University Taught by: Jeff Leek, PhD

## Week 1

1. Reproducibe research

Case study:
["Genomic signatures to guide the use of chemotherapeutics"](https://www.nature.com/articles/nm1491)
-> claim that you could use genomic signatures (gene expression signatures) to guide the use of chemotherapeutics. So that you could tell which chemotherapeutic would work best for which person based on their gene expression profiles.

However,

- many of data not avaiable

- many of code not avaiable

- once reserchers got the data and code, it turned out it did not produce the results they had in the paper

- if you runned the code on different days, times, you got different answers (there was random component to the predicitons)

But, some clinical trials had been launched based on paper. People were assigned to chemotherapy possibly erroneously based on original analysis. Led to law suits.

2.  What you need to share research anlysis

    - The raw data set

    - A tidy data set

    - A code book describing each variable and its values in the tidy dataset (technology, units, machines...)

    - An explicit and exact recipe you used to go from 1 -> 2,3 (Code that get raw data and otputs tidy one either R or more preferably R markdown or iPython notebook)

3.  Good practise

    Always add:

            ```{r session info}
            sessionInfo()
            ```

    &

            The document was processed on `r Sys.Date()`.

    at the end of your notebooks.

4.  The three table in genomics

<figure>
  <img src="./images/summarized.png" alt="From [1]">
  <figcaption>Figure 1: From [1]</figcaption>
</figure>

5. Exploratory analysis

## Week 2

### Visualizing data

- often thousands features per sample, so to visualize/see patterns dimetion reduction is nice

  - SVD genes x samples decomposed into:

    - U eigenarrays/left sigular vectors: patterns accross the arrays
    - D singular values - diagonal matrix (signifies how much of the variance is explained by those various patterns ^)
    - V transposed: eigengenes: relationship in column patterns - so patterns across genes

    Properties:

    - columns od Vt/rows of U are orthogonal (uncorellated with each other)

  - Other: Multidimentional scaling, independent components analysis, non-negative matrix factorization
  - More on the topic [here](https://courses.edx.org/courses/course-v1:HarvardX+PH525.3x+1T2018/0b42cffa7c6e4c559bf74f93fb864a59/).
  - prcomp: Principal Component Analysis

- svd, boxplots, outlier analysis

### Making data comparable (Pre-processing & normalization)

- Quantile normalization: force the distributions to be exactly the same

<figure>
  <img src="./images/quantile.png" alt="From [1]">
  <figcaption>Figure 2: From [2]</figcaption>
</figure>

- Preprocessing and normaliztion are **higly platform/problem dependent** :(

### Linear modelling

Full course on that [here](https://www.coursera.org/learn/regression-models).

- case: a continuous outcome, but maybe a not continuous covariate or a categorical covariate or a factor-level covariate

- Adjusting for covariates in linear regression models.

### Batch-effects, Confounders

### Citations

[1] Huber, W., Carey, V., Gentleman, R. et al. Orchestrating high-throughput genomic analysis with Bioconductor. Nat Methods 12, 115–121 (2015). https://doi.org/10.1038/nmeth.3252

[2] Hicks, S.C., Irizarry, R.A. quantro: a data-driven approach to guide the choice of an appropriate normalization method. Genome Biol 16, 117 (2015). https://doi.org/10.1186/s13059-015-0679-0

### To read

1. [Points of View columns on data visualization published in Nature Methods](https://blogs.nature.com/methagora/2013/07/data-visualization-points-of-view.html)

2. [Using Microsoft Excel to obscure your data and annoy your readers](https://www.biostat.wisc.edu/~kbroman/presentations/IowaState2013/graphs_combined.pdf)
