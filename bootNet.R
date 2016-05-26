#!/usr/bin/Rscript
#
# Created: 2016/01/10
# Last modified: 2016/03/03
# Author: Miles Benton
# Version: 0.1.1.1
# 
# bootNet is a wrapper for the fantastic glmnet R package - it brings bootstrapping and parallel processing to the elastic-net framework.
#
# For more information please see README at: https://github.com/sirselim/bootNet
#
# """
#
# This update adds a parallel version of the bootNet function to utalise multiple cores if available
# WARNING [here be dragons!]: be aware of the amount of available system RAM when using bootNet.parallel, if the data
# set is large even running across 4-8 cores will quickly utalise many GB of RAM - you have been warned!
#
# This script is currently set up to analyse methylation data in the form of beta matrices.
# The beta matrix must have CpG sites as rows and samples as columns for bootNet to work.
# Hopefully future updates will allow analysis of various types of data.
#
# E.g. use:
# bootNet.parallel(data = x, outcome = y, Alpha = 0.1, iter = 1000, sub_sample = 0.666, cores = 4, sampleID = sampleID)
#
# """

########################
## bootstrap function ##
########################
bootNet <- function(data, outcome, Alpha, iter, Lambda, sub_sample, sampleID){
  
  # report on outcome type
  if (is.numeric(outcome) == TRUE) {
    cat('...outcome is quantitative, using gaussian approach in glmnet model...', '\n')
  } else if (is.factor(outcome) == TRUE) {
    cat('...outcome is qualitative, using binomial approach in glmnet model...', '\n')
  } else {
    outcome.error <- paste0('outcome is [', class(outcome), ']', ' it needs to be numeric or factor...')
    cat('...ERROR:', outcome.error, '\n')
  }
  
  ## implement a check for NA's in data and outcome, quit with error if found
  outcome_nas <- any(is.na(outcome))
  data_nas <- any(is.na(data))
  if (!outcome_nas | !data_nas) stop("bootNet has discovered NA's in outcome or data, terminating function call") 

  # load packages
  require(glmnet)
  
  # create empty list
  cpg_list <- list()
  
  # create an object to extract sites from later
  cpg.sites <- rownames(data)
  
  # transpose data for glmnet
  data <- t(data)
  
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
      # Do glmnet - 'gaussian' family
      fit <- glmnet(x = newData, y = newOut, family = "gaussian", alpha = Alpha)
    } else {
      # sample from qualitative outcome
      newDataInd <- c(sample(grep(levels(outcome)[1], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[1], outcome))))), 
                      sample(grep(levels(outcome)[2], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[2], outcome))))))
      newData <- data[newDataInd,]  # subset the data
      newOut <- outcome[newDataInd] # In the outcome variable get the same patients as were selected for this iteration
      # Do glmnet - 'binomial' family
      fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    } 
    
    # Get model coefficients
    Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
    
    # Get CpG list for which coefficients are not 0
    cpgs <- cpg.sites[Coefficients@i]
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
bootNet.parallel <- function(data, outcome, Alpha, iter, Lambda, sub_sample, cores, sampleID){

  # report on outcome type
  if (is.numeric(outcome) == TRUE) {
    cat('...outcome is quantitative, using gaussian approach in glmnet model...', '\n')
  } else if (is.factor(outcome) == TRUE) {
    cat('...outcome is qualitative, using binomial approach in glmnet model...', '\n')
  } else {
    outcome.error <- paste0('outcome is [', class(outcome), ']', ' it needs to be numeric or factor...')
    cat('...ERROR:', outcome.error, '\n')
  }

  ## implement a check for NA's in data and outcome, quit with error if found
  outcome_nas <- any(is.na(outcome))
  data_nas <- any(is.na(data))
  if (!outcome_nas | !data_nas) stop("bootNet has discovered NA's in outcome or data, terminating function call") 
  
  # load packages
  require(glmnet)
  require(foreach)
  require(doParallel)
  
  # register cores
  registerDoParallel(cores = cores)
  
  # create an object to extract sites from later
  cpg.sites <- rownames(data)
  
  # transpose data for glmnet
  data <- t(data)
  
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
      # Do glmnet - 'gaussian' family
      fit <- glmnet(x = newData, y = newOut, family = "gaussian", alpha = Alpha)
    } else {
      # sample from qualitative outcome
      newDataInd <- c(sample(grep(levels(outcome)[1], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[1], outcome))))), 
                      sample(grep(levels(outcome)[2], outcome), ceiling(sub_sample*(length(grep(levels(outcome)[2], outcome))))))
      newData <- data[newDataInd,]  # subset the data
      newOut <- outcome[newDataInd] # In the outcome variable get the same patients as were selected for this iteration
      # Do glmnet - 'gaussian' family
      fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    } 
    
    # Get model coefficients
    Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
    # Get CpG list for which coefficients are not 0
    cpgs <- cpg.sites[Coefficients@i]
    list(cpgs)
  }
}
#################################
