#!/usr/bin/Rscript
#
# Created: 2016/01/10
# Last modified: 2016/03/02
# Author: Miles Benton
# Version: 0.1.1.0
# 
# This update add a parallel version of the bootNet function to utalise multiple cores if available
# WARNING: be aware of the amount of available system RAM when using bootNet.parallel, if the data
# set is large even running across 4-8 cores will quickly utalise many GB of RAM - you have been warned!
#
# """
# This script pulls out information from vcf files.
# The script accepts 1 argument (file name).
#
# E.g. INPUT
# ./vcfcompiler_diagnostics.sh /path/to/initial_filter.vcf
#
# The output is a csv file named after the vcf
# WARNING: if the csv file already exists, this script will REPLACE ALL OF ITS CONTENT
# E.g. OUTPUT
# initial_filter.csv
# """

########################
## bootstrap function ##
########################
bootNet <- function(data, outcome, Alpha, iter, Lambda, beta_matrix, sub_sample){
  # load packages
  require(glmnet)
  # require(survival)
  cpg_list <- list()
  # bootstrap process 
  for (i in 1:iter){
    set.seed(i)
    # Select a random number of patients
    # this is currently only set up for 30 samples (15 control vs 15 case), need to fix this
    # newDataInd <- c(sample(1:15,floor(sub_sample*15)), sample(16:30,floor(sub_sample*15))) # sampling from both groups (quantitative)
    newDataInd <- c(sample(grep('Blood', outcome), floor(sub_sample*(length(grep('Blood', outcome))/2))), 
                    sample(grep('Buccal', outcome), floor(sub_sample*(length(grep('Buccal', outcome))/2))))
    # above is hard coded - need to implement a check for outcome (y) type, i.e. qualitative vs quantitative
    # then perform the correct family model
    #Subset the data
    newData <- data[newDataInd,]
    # In the outcome variable get the same patients as were selected for this iteration
    newOut <- outcome[newDataInd]
    # Do glmnet
    fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    # Get model coefficients for glmnet
    Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
    # Get CpG list for which coefficients are not 0
    cpgs <- rownames(beta_matrix[Coefficients@i,])
    name <- paste('run:', i, sep = '')
    cpg_list[[name]] <- cpgs
    print(i)
  }
  return(cpg_list)
  cat('\n', ' ...Processing Done...')
}
########################

#################################
## parallel bootstrap function ##
#################################
bootNet.par <- function(data, outcome, Alpha, iter, Lambda, beta_matrix, sub_sample, cores){
  # load packages
  require(glmnet)
  library(foreach)
  library(doParallel)
  # register cores
  registerDoParallel(cores = cores)
  # bootstrap process 
  foreach (i = 1:iter, .combine = c) %dopar% {
    set.seed(i)
    # Select a random number of patients/samples
    newDataInd <- c(sample(grep('Blood', outcome), floor(sub_sample*(length(grep('Blood', outcome))/2))), 
                    sample(grep('Buccal', outcome), floor(sub_sample*(length(grep('Buccal', outcome))/2))))
    #Subset the data
    newData <- data[newDataInd,]
    # In the outcome variable get the same patients as were selected for this iteration
    newOut <- outcome[newDataInd]
    # Do glmnet
    fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    # Get model coefficients for glmnet
    Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
    # Get CpG list for which coefficients are not 0
    cpgs <- rownames(beta_matrix[Coefficients@i,])
    list(cpgs)
  }
}
#################################