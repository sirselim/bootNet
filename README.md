# bootNet

`bootNet` is a wrapper for the fantastic [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) `R` package - it brings bootstrapping and parallel processing to the elastic-net framework.

This script is currently set up to analyse methylation data in the form of beta matrices. The beta matrix must have CpG sites as rows and samples as columns for bootNet to work. Hopefully future updates will allow analysis of various types of data.

## Updates

### Version: 0.1.1.1

This update adds a parallel version of the bootNet function to utalise multiple cores if available.

**WARNING:** be aware of the amount of available system RAM when using `bootNet.parallel()`, if the data set is large even running across 4-8 cores will quickly utalise many GB of RAM - **you have been warned!**

## Example usage:

What do you need to provide the `bootNet` script?

  - *data* - methylation beta matrix (formated CpG sites as rows, samples as columns)
  - *outcome* - either as numeric (quantitative trait) or as factor (qualitative trait)
  - *Alpha* - sets the mixture between ridge-regression and lasso
  - *iter* - the number of iterations to bootstrap
  - *sub_sample* - percentage of sample/case-control groups to sub-sample for bootstrapping
  - *sampleID* - a list of sample IDs (usually the column names of the beta matrix). These **MUST** be in the same order as the samples in the beta matrix
  - *cores* - number of cores to use for `bootNet.parallel`

### bootNet()
`bootNet(data = x, outcome = y, Alpha = 0.1, iter = 1000, sub_sample = 0.666, sampleID = sampleID)`

### bootNet.parallel()
`bootNet.parallel(data = x, outcome = y, Alpha = 0.1, iter = 1000, sub_sample = 0.666, cores = 4, sampleID = sampleID)`

## To do list

  - add example data set
  - add a decent tutorial
  - include additional quality checks
  - explore bias estimation and performance of selected markers
  - implement a double-bootstrap method as an additional way to explore performance metrics
    + separate training and test data
    + bootstrap both data sets
    + compare selected sites
  - look into AUC and ROC as another form of marker selection/validation  

## Performance expectations

A few test examples showing systems and performance metrics.

**Experiment 1**: 24 samples run on Illumina 450K methylation array (24 columns, 446280 rows). Phenotype was quantitative (age).

System used:

  - Intel Core i7-4900MQ (4 core, 8 threads)
  - 32GB DDR3
  - 256GB SSD
  - Linux - 4.3.0-towo.3-siduction-amd64 x86_64 (64 bit)

Running **100 iterations** and using **4 cores** the maximum observed RAM usage was **20GB**, taking **85.6 seconds**. 

**Experiment 2**: 75 samples run on Illumina 450K methylation array (24 columns, 445998 rows). Phenotype was quantitative (age).

  - Intel Xeon E5-2620 (6 cores, 12 threads)
  - 128GB DDR3
  - 256GB SSD
  - Linux - Debian Sid (64 bit)

Running **100 iterations** and using **6 cores** the maximum observed RAM usage was **24GB**, taking **145 seconds**. 

## Dependencies 

There are 3 R packages required by `bootNet`:

  - `glmnet`
  - `foreach` (bootNet.parallel)
  - `doParallel` (bootNet.parallel)