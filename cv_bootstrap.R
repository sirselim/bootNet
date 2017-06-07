####
## bootstrapped cross-validation
## Miles Benton
## date created: 1st Nov 2016
## date modified: 21st Apr 2017 

## this is an experimental script to explore generated bootstrapped glmnet CV runs
## once the set number of iterations has been performed a matrix of 2 columns is returned
## these columns represent the lambda and 1SE lambda from each CV run

##
# this only needs to be run once on a data set (takes a little while)
# BEWARE: if running in parallel large RAM usage (>26GB for this data)
require(glmnet)
require(doMC)
registerDoMC(cores=4)
# NOTE: cross val doesn't like factors for y...

# picked a random (30) number of iterations (CV performs a certain number itself, which is why this takes a little time)
lambda.matrix <- NULL
for (i in c(1:30)) {
  
  # if y is quantitative use gaussian
  # cv.fit <- cv.glmnet(x, as.numeric(y), parallel = T)  # Perform cross validation on the fited model
  # if y is qualitative use binomial
  cv.fit <- cv.glmnet(t(snp.matrix), y, family = 'binomial', parallel = T)  # Perform cross validation on the fited model
  # create data.frame of lambda.min and lambda.1se
  lambda.matrix <- rbind(lambda.matrix, data.frame(lambda.min = cv.fit$lambda.min, lambda.1se = cv.fit$lambda.1se))
  
}

# take the mean lambdas
colMeans(lambda.matrix)	# user can decide which 'mean' lambda to use moving into bootnet
##