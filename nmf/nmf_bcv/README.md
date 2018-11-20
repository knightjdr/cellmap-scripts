# nmf-bcv: a workflow for non-negative matrix factorization (NMF) bi-cross-validation (BCV)

These Python/R scripts are meant to accompany a manuscript that will be submitted in future. This workflow implements bi-cross-validation with simple residuals from the following paper: Owen, A. B. & Perry, P. O. [Bi-cross-validation of the SVD and the nonnegative matrix factorization](https://arxiv.org/abs/0908.2062). *The annals of applied statistics* (2009).

**Important: As the manuscript for this code has not yet been published, please do not redistribute it without prior permission from [Simon Eng](mailto:simon.eng@mail.utoronto.ca).**

## Dependencies

-   Python 3.5 or higher, with the following packages installed:
    +   feather (`pip install feather-format`)
    +   numpy
    +   pandas
    +   scikit-learn (`pip install scikit-learn`)
    +   scipy
    +   tqdm (`pip install tqdm`)

## Input file

A comma-separated values (CSV) file with observations, patients, etc. as rows and measurements as columns. The first column should be an identifier for each row.

We recommend scaling the input data to unit variance. (JK note: I don't think this should be done for spectral count data)

## Workflow

### 1. BCV

        python cv_nmf_bcv.py --input input.csv --output q2.feather

Options:

`--init`
: Method to initialize values with. Possible values include `random`, `nndsvd` (non-negative double singular value decomposition), `nndsvda`, and `nndsvdar`

`--l1-ratio`
: Regularization mixing parameter. Use `0` for an L2 penalty, `1` for an L1 penalty, or any number between 0--1 for a combination

`--k`
: Factorization ranks to test (e.g., `--k 1 2 3 4 5 6 7 8 9 10`)

`--folds`
: The number of folds to split the input data into.

`--iterations`
: The number of iterations to run bi-cross-validation with

`--seed`
: The seed to use to initialize bi-cross-validation

`--cores`
: The number of processor cores to use *(you will want to specify this in most cases)*

### 2. Obtaining the optimal factorization rank

        python get_bcv_rank.py --input q2.feather --output k.txt

Options:

`--q2-threshold`
: Some runs of BCV may produce Q2 << -1, which will skew results. By default, Q2 values less than -1 will be omitted from analysis

### 3. NMF

        python nmf.py --input input.csv --k `cat k.txt` --model-output model.pkl --basis-output basis.csv --score-output scores.csv

*Note: Scores are often also referred to as "coefficients".*

Options:

`--init`
: See the [BCV](#1-bcv) section.

`--l1-ratio`
: See the [BCV](#1-bcv) section.

### 4. Obtaining cluster assignments

        python get_clusters.py --input scores.csv --output clusters.csv

Options:

`--allow-unassigned`
: Permits observations to not be assigned to any cluster

`--variance-coefficient`
: Requires coefficients to be above a given number of standard deviations above 0 to receive a cluster assignment. This option exists primarily for observations that cannot be explained by any factors
