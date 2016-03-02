#!/usr/bin/Rscript
#
# Created: 2016/01/10
# Last modified: 2016/03/02
# Author: Miles Benton
# Version: 0.1.1.1
# 
# This update add a parallel version of the bootNet function to utalise multiple cores if available
# WARNING: be aware of the amount of available system RAM when using bootNet.parallel, if the data
# set is large even running across 4-8 cores will quickly utalise many GB of RAM - you have been warned!
#
# """
#
# This script pulls out information from vcf files.
# 
# E.g. use:
# bootNet.parallel(data = x, outcome = y, Alpha = 0.1, iter = 1000, beta_matrix = beta_norm, sub_sample = 0.666, cores = 4, sampleID = sampleID)
#
# """

########################
## bootstrap function ##
########################
bootNet <- function(data, outcome, Alpha, iter, Lambda, beta_matrix, sub_sample, sampleID){
  # include checks for type of outcome data
  # quantitative needs to be named numeric
  # qualitative needs to be factor
  
  # load packages
  require(glmnet)
  
  # create empty list
  cpg_list <- list()
  
  # bootstrap process 
  for (i in 1:iter){
    set.seed(i)
    # Select a random sub-sample from all samples
    # first determine whether outcome is qualitative or quantitative
    if (is.numeric(outcome) == TRUE) {
      # sample from quantitative outcome
      # NOTE: sampleID must be the same order as samples in the beta matrix/data!
      subID <- sample(sampleID, ceiling(sub_sample*(length(sampleID))))  # get ID's for sub-sample
      newDataInd <- outcome[names(outcome) %in% subID]  # subset outcome for correct samples
      newData <- data[rownames(data) %in% names(newDataInd),] # subset the data
      newOut <- as.numeric(newDataInd)
    } else {
      # sample from qualitative outcome
      newDataInd <- c(sample(grep(levels(outcome)[1], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[1], outcome))))), 
                      sample(grep(levels(outcome)[2], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[2], outcome))))))
      newData <- data[newDataInd,]  # subset the data
      newOut <- outcome[newDataInd] # In the outcome variable get the same patients as were selected for this iteration
    } 
    
    # Do glmnet
    # determine whether outcome is qualitative or quantitative
    if (is.numeric(outcome) == TRUE) {
      # if TRUE then 'gaussian' family
      cat('...outcome is quantitative, using gaussian approach in glmnet model...')
      fit <- glmnet(x = newData, y = newOut, family = "gaussian", alpha = Alpha)
    } else {
      # if FASLE then 'binomial' family
      cat('...outcome is qualitative, using binomial approach in glmnet model...')
      fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    }
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
bootNet.par <- function(data, outcome, Alpha, iter, Lambda, beta_matrix, sub_sample, cores, sampleID){
  # load packages
  require(glmnet)
  library(foreach)
  library(doParallel)
  # register cores
  registerDoParallel(cores = cores)
  # bootstrap process 
  foreach (i = 1:iter, .combine = c) %dopar% {
    set.seed(i)
    
    # Select a random sub-sample from all samples
    # first determine whether outcome is qualitative or quantitative
    if (is.numeric(outcome) == TRUE) {
      # sample from quantitative outcome
      # NOTE: sampleID must be the same order as samples in the beta matrix/data!
      subID <- sample(sampleID, ceiling(sub_sample*(length(sampleID))))  # get ID's for sub-sample
      newDataInd <- outcome[names(outcome) %in% subID]  # subset outcome for correct samples
      newData <- data[rownames(data) %in% names(newDataInd),] # subset the data
      newOut <- as.numeric(newDataInd)
    } else {
      # sample from qualitative outcome
      newDataInd <- c(sample(grep(levels(outcome)[1], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[1], outcome))))), 
                      sample(grep(levels(outcome)[2], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[2], outcome))))))
      newData <- data[newDataInd,]  # subset the data
      newOut <- outcome[newDataInd] # In the outcome variable get the same patients as were selected for this iteration
    } 
    
    # Do glmnet
    # determine whether outcome is qualitative or quantitative
    if (is.numeric(outcome) == TRUE) {
      # if TRUE then 'gaussian' family
      cat('...outcome is quantitative, using gaussian approach in glmnet model...')
      fit <- glmnet(x = newData, y = newOut, family = "gaussian", alpha = Alpha)
    } else {
      # if FASLE then 'binomial' family
      cat('...outcome is qualitative, using binomial approach in glmnet model...')
      fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    }
    
    # Get model coefficients for glmnet
    Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
    # Get CpG list for which coefficients are not 0
    cpgs <- rownames(beta_matrix[Coefficients@i,])
    list(cpgs)
  }
}
#################################