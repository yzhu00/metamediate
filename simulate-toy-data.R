## Simulate toy data to use as an example 
## Effect sizes: large effects from simulations 
## 10 studies, as a list
## Dataframes following the same format: A (treatment), M1, M2, covar1, covar2, y (outcome)

generateToyData <- function(randomseed){

  
  ## effects 
  
  # realbeta.m1 <- c(0.03, 0.1, 0.004, -0.004) ## trt on weight gain = 0.1 
  realbeta.m1 <- c(3, 3, 4, -4)
  names(realbeta.m1) <- c("intercept", "A", "cov1", "cov2")
  
  realbeta.m2 <- c(7, -5, 0.5, 0.1)
  names(realbeta.m2) <- c("intercept", "A", "cov1", "cov2")
  
  ## positive interaction between weight gain and treatment
  realbeta.y <- c(5, -2, -0.5, 0.6, -1, 1, 0.5, -0.5, 0.01, 0.03) 
  
  names(realbeta.y) <- c("intercept", "A", "M1", "M2",
                         "cov1", "cov2", "A:M1", "A:M2", 
                         "M1:M2", "A:M1:M2")
  
  
  #### Simulate data using the set parameters ####
  
  ## two covariates, bivariate normal (basically PANSS+ baseline, PANSS- baseline)
  ## numbers are observed from the merged real data before imputation 
  covars <- mvrnorm(popN, mu = c(24, 23.5), Sigma = matrix(c(23, -3, -3, 30), nrow =2))
  
  ## treatment: Bernoulli p = 0.5 
  A <- rbinom(popN, size = 1, prob = 0.5)
  
  ## mediator 1: weight gain 
  ## changing the definition to percent weight gain
  mu1 <- cbind(rep(1, popN), A, covars) %*% realbeta.m1
  M1 <- rnorm(popN, mean = mu1, sd = sqrt(0.5)) # larger variance than the observed data 
  
  #pgain <- mu + rnorm(n = popN, mean = 0, sd = sqrt(0.001))
  
  ## mediator 2: PANSS + at week 4/5 
  mu2 <- cbind(rep(1, popN), A, covars) %*% realbeta.m2
  ## PANSS scores are technically bounded by [7,49]
  ## simulate as a truncated normal variable 
  M2 <- rnorm(popN, mean = mu2, sd = sqrt(24))
  
  
  ## outcome: PANSS - at week 6 
  mu3 <- cbind(rep(1, popN), A, M1, M2, covars, 
              A*M1, A*M2, M1*M2, A*M1*M2) %*% realbeta.y 
  y <- rnorm(popN, mean = mu3, sd = sqrt(10))
  
  ## assemble the populatuion dataset 
  population <- as.data.frame(cbind(A, M1, M2, covars, y))
  colnames(population) <- c("A", "M1", "M2", "covar1", "covar2", "y")## Compute the effects in the population 
  
  
  study.data <- lapply(1:nstudy, FUN = function(x){
    samplesize <- 200 # keeping it realistic
    study <- population[sample(1:popN, size = samplesize, replace = FALSE),]
  }
  )
  
  return(study.data)
}



