# Created: 2020/07/22
# Last modified: 2020/07/22
# Author: Sean Burnard
# Version: 0.2.1
#
# Creates a function on top of BootNet (a wrapper for GLMNet developed by Miles Benton) 
#   for SNP analysis from recoded plink files
# This extracts all SNP selected during the BootNet procedure (iterations of GLMNet)
#   And converts them into percentages while allowing (and recording) the optional variables
# Iterations of the cv.fit procedure (to obtain lambda) is also advised (with 30 recommended and set by default). However, this can take a long time with larger dataset. 
#   If you do not wish to perform iterations of cv.fit this can simply be set to 1.
#   Iterations of cv.fit was introduced as there is some variability each time cv.fit is run due to the way the data is split during the glmnet cross validaiton procedure. 
#   Therefore, this provides some stabilisation around a consistent (average) lambda.
#
## Raw file obtained by plink using --recode A
## CSV file name will automatically append ".csv" (no need to write .csv at the end)
## Default GLMNet iterations is 200, this can be changed as desired. Advise putting the number used within the file name!
## Start with a low number of BootNet iterations then increase (beyond the number of samples) and until the 'top hits' appear to settle i.e. if variables are being selected by GLMNet above 5% then the highest variables tend to stabalise in order as the number of iterations is increased and this may also start to tease out any set of variables that may initially be selected at similar freqeuncies (especially those reaching ~100%).



SNPs_For_BootNet <- function(RAW_FILE, CSV_FILE_NAME,
                               Number_of_GLMNet_Iterations = 200, CV_Iterations = 30,
                               a_alpha = 0.1, Subsample = 0.66,
                               penalty = "min") ## penalty can be either "min" or "1se" (referring to lambda.min or lambda.1se, respectively)
  { 

# Load dataset and extract SNP matrix and phenotypes with IDs (in order).

  AllCases_Recoded <- read.table(RAW_FILE, header = T)
  AllCases_Recoded2 <- AllCases_Recoded
  AllCases_Recoded2[1:6] <- list(NULL)
  AllCases_Recoded2 <- as.matrix(AllCases_Recoded2) #convert to matrix
  AllCases_IDs <- AllCases_Recoded$IID #Obtain IDs
  row.names(AllCases_Recoded2) <- AllCases_IDs #Add ID names in order to SNP matrix


  AllCases_Pheno <- as.factor(AllCases_Recoded$PHENOTYPE) #Obtain pheno (in same order) and make factor!!

###Load GLMNet and other required packages
  library("glmnet")
  library(dplyr)
  
### Bootstrapp of cvfit to use for GLMNet BootNet run

  # picked a random (30) number of iterations (CV performs a certain number itself, which is why this takes a little time)
  lambda.matrix <- NULL
 
   for (i in c(1:CV_Iterations)) {
  
  # if y is quantitative use gaussian
  # cv.fit <- cv.glmnet(x, as.numeric(y), parallel = T)  # Perform cross validation on the fited model
  # if y is qualitative use binomial
    cv.fit <- cv.glmnet(AllCases_Recoded2, AllCases_Pheno, family = "binomial", alpha = a_alpha)  # Perform cross validation on the fited model
  # create data.frame of lambda.min and lambda.1se
    lambda.matrix <- rbind(lambda.matrix, data.frame(lambda.min = cv.fit$lambda.min, lambda.1se = cv.fit$lambda.1se))
  
  }

# take the mean lambdas
  print(colMeans(lambda.matrix))	# user can decide which 'mean' lambda to use moving into bootnet

  lamda.Means <- colMeans(lambda.matrix)


### Run GLMNet with BootNet script from Miles  
  
  source("./bootNet.R")
  BootNet_RUN <- bootNet(
    data = t(AllCases_Recoded2), #  transposes matrix so col = samples, rows = SNPs
    outcome = AllCases_Pheno, 
    Alpha = a_alpha,
    iter=Number_of_GLMNet_Iterations,
    sub_sample = Subsample,
    sampleID = AllCases_IDs, 
    Lambda = if(penalty == "min") {
                  colMeans(lambda.matrix)[["lambda.min"]]
                  }
              else if (penalty == "1se") {
                  colMeans(lambda.matrix)[["lambda.1se"]]
              }
              else {
                print("missing correct lambda choice")
              }
      
    ) ## if no bootstrapp for cv  -->  cvfit_AllCases$lambda.min
  ## method = "BOOTSTRAP",    ##No method option for parrallel

  
  
## Obtain some meaningful results from the BootNet run and create csv file for results and text file for settings
  BootNet_RUN_table <- data.frame(snp = names(sort(table(unlist(BootNet_RUN)))),
                          iterations = as.numeric(sort(table(unlist(BootNet_RUN)))))

  BootNet_RUN_table <- BootNet_RUN_table[order(-BootNet_RUN_table$iterations),]
  
  BootNet_RUN_table <- BootNet_RUN_table %>% mutate(Perc_IT = (iterations / Number_of_GLMNet_Iterations * 100))
  

  write.csv(BootNet_RUN_table, file = paste0(CSV_FILE_NAME, ".csv"))
  
  
  
  Settings_used <- data.frame(
                    Variable = c("Dataset:",
                             "alpha:",
                             "Subsample:",
                             "CV_iterations:",
                             "GLMNet_iterations:",
                             "lambda_used:",
                             "lambda.min:",
                             "lambda.1se:"
                             ), 
                    Value = c(RAW_FILE,
                          a_alpha,
                          Subsample,
                          CV_Iterations,
                          Number_of_GLMNet_Iterations,
                          penalty,
                          colMeans(lambda.matrix)[["lambda.min"]],
                          colMeans(lambda.matrix)[["lambda.1se"]]
                          ))
    
  write.table(Settings_used, file = paste0(CSV_FILE_NAME, "_Settings_Used.txt"), row.names = F, quote = F)
  
  
  paste0("written values used to '", CSV_FILE_NAME, ".txt'",
        "; and results written to '", CSV_FILE_NAME, ".csv' " )
  
}