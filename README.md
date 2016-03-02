# bootNet

`bootNet` is a wrapper for the fantastic [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) `R` package - it brings bootstrapping and parallel processing to the elastic-net framework.

This script is currently set up to analyse methylation data in the form of beta matrices. Hopefully future updates will allow analysis of various types of data.

# Updates

## Version: 0.1.1.1

This update adds a parallel version of the bootNet function to utalise multiple cores if available

**WARNING:** be aware of the amount of available system RAM when using `bootNet.parallel()`, if the data set is large even running across 4-8 cores will quickly utalise many GB of RAM - **you have been warned!**

# Example usage:

## bootNet()
`bootNet.parallel(data = x, outcome = y, Alpha = 0.1, iter = 1000, sub_sample = 0.666, cores = 4, sampleID = sampleID)`

## bootNet.parallel()
`bootNet.parallel(data = x, outcome = y, Alpha = 0.1, iter = 1000, sub_sample = 0.666, cores = 4, sampleID = sampleID)`

## To do list

  - add example data set
  - add a decent tutorial
  - include additional quality checks

