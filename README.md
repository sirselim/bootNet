# bootNet ('strapping the elastic-net)

`bootNet` is a wrapper for the fantastic [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) `R` package - it brings bootstrapping and parallel processing to the elastic-net framework.

This script was originally desinged to analyse methylation data in the form of beta matrices. The beta matrix must have CpG sites as rows and samples as columns for bootNet to work. 

With a recent update `bootNet` is now data agnostic and has been tested with SNP and expression data (as well as methylation).

## Updates

### Version: 0.1.2.0

  - code base overhauled to be data agnostic (i.e. it can be used for methylation, SNP, expression data, but also any matrix of data if set up correctly in normal `glmnet` format)
  - added an example data set in `/data`, this is a collection of blood and buccal cell samples 

### Version: 0.1.1.1

  - This update adds a parallel version of the bootNet function to utilise multiple cores if available.

**WARNING:** be aware of the amount of available system RAM when using `bootNet.parallel()`, if the data set is large even running across 4-8 cores will quickly utalise many GB of RAM - **you have been warned!**  Some real-world usage metrics across different Linux systems are provided below.

  - Added functionality to detect NAs in data/outcome and terminate the function

  - Created a gh-pages branch that will eventually house a website for `bootNet`: [here](https://sirselim.github.io/bootNet/)

## Example usage:

What do you need to provide the `bootNet` script?

  - *data* - data matrix (formated predictors as rows, samples as columns)
  - *outcome* - either as numeric (quantitative trait) or as factor (qualitative trait)
  - *Alpha* - sets the mixture between ridge-regression and lasso
  - *iter* - the number of iterations to bootstrap
  - *sub_sample* - percentage of sample/case-control groups to sub-sample for bootstrapping
  - *sampleID* - a list of sample IDs (usually the column names of the beta matrix). These **MUST** be in the same order as the samples in the beta matrix
  - *method* - optional string specifying one of two methods ('JACKKNIFE', 'LOOCV')
  - *cores* - number of cores to use for `bootNet.parallel`

### bootNet()
`bootNet(data = x, outcome = y, Alpha = 0.1, iter = 1000, sub_sample = 0.666, sampleID = sampleID, method = method)`

### bootNet.parallel()
`bootNet.parallel(data = x, outcome = y, Alpha = 0.1, iter = 1000, sub_sample = 0.666, cores = 4, sampleID = sampleID)`

## To do list

  - add functionality to automatically 'tally' the results and rank them by number of iterations passed 
  - ~~merge development script back into the main code base~~
  - ~~add BOOTNET, JACKKNIFE and LOOCV methods to main script~~
  - add BOOTNET, JACKKNIFE and LOOCV methods to the parallel method
  - add example data set and a decent tutorial
  - include additional quality checks
    + ~~ensure outcome is either numeric or factor~~
    + ~~check for missing data (NAs) in both data and outcome~~
  - explore bias estimation and performance of selected markers
  - implement a double-bootstrap method as an additional way to explore performance metrics
    + separate training and test data
    + bootstrap both data sets
    + compare selected sites
  - look into AUC and ROC as another form of marker selection/validation
  - ~~add an argument to perform "leave one out" analysis (alternative to sub-sampling)~~
  - write up as manuscript
  - add option for paired samplesin LOO-CV method
    + ***bootNet: identifying robust classifiers in methylation data***  

## Performance expectations

A few test examples showing systems and performance metrics.

**Experiment 1**: 24 samples run on Illumina 450K methylation array (24 columns, 446280 rows). Phenotype was quantitative (age).

System used (Dell laptop: Precision M4800):

  - Intel Core i7-4900MQ (4 core, 8 threads)
  - 32GB DDR3
  - 256GB SSD
  - Linux - 4.3.0-towo.3-siduction-amd64 x86_64 (64 bit)
    + Running **100 iterations** and using **4 cores** the maximum observed RAM usage was **20GB**, taking **85.6 seconds**. 

**Experiment 2**: 75 samples run on Illumina 450K methylation array (75 columns, 445998 rows). Phenotype was quantitative (glucose).

System used (Dell workstation):

  - Intel Xeon E5-2620 (6 cores, 12 threads)
  - 128GB DDR3
  - 256GB SSD
  - Linux - Debian Sid (64 bit)
    + Running **100 iterations** and using **10 cores** the maximum observed RAM usage was **30GB**, taking **120.9 seconds**.  
    + Running the above at **1000 iterations** and using **10 cores** the maximum observed RAM usage was **30GB**, taking **1152.8 seconds**. 
    + Running the above at **5000 iterations** and using **10 cores** the maximum observed RAM usage was **30GB**, taking **5718.8 seconds**. 

**Experiment 3**: 105 samples run on Illumina 450k methylation array (105 columns, 380777 rows). Phenotype was qualitative (BMI separated into two categories)

System used (Z170-D3H)

  - Intel Core i7-6700 3.40GH (4 cores, 8 threads)
  - 16GB DDR4 (2 x 8GB @ 2133 MHz)
  - 2TB Western Digital Harddrive
  - Linux Ubuntu 14.04 LTS (64 bit)
    + **10 iterations** and using **1 core** the maximum observed RAM usage was **5GB**, taking **34.42089 secs**.
    + **10 iterations** and using **2 cores** the maximum observed RAM usage was **9GB**, taking **20.0074 secs**.
    + **10 iterations** and using **3 cores** the maximum observed RAM usage was **12.5GB**, taking **18.31608 secs**.
    + **20 iterations** and using **3 cores** the maximum observed RAM usage was **13.8GB**, taking **33.594 secs**.
    + **30 iterations** and using **3 cores** the maximum observed RAM usage was **12.4GB**, taking **48.210 secs**.
    + **50 iterations** and using **3 cores** the maximum observed RAM usage was **13GB**, taking **83.090 secs**.
    + **100 iterations** and using **3 cores** the maximum observed RAM usage was **13.6GB**, taking **161.983 secs**.
    + **1000 iterations** and using **3 cores** the maximum observed RAM usage was **13 GB**, taking **1628.748 secs**.


## Dependencies 

There are 3 R packages required by `bootNet`:

  - `glmnet`
  - `foreach` (bootNet.parallel)
  - `doParallel` (bootNet.parallel)
