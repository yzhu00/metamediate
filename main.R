
## main script 
## showing how to run a meta-analysis of decomposed causal effects 
## using simulated toy data

## load packages 
library(MASS)
library(mvmeta)
library(Matrix)
library(msm)


# Part 1: simulate toy data -----------------------------------------------

## homogenous effects across trials 

## set parameters 
nstudy <- 10 
popN <- 1000000 # population size 

source("simulate-toy-data.R")
studies <- generateToyData(randomseed = 42) 

# Part 2: perform regression analysis in each study -----------------------

## fit regression models in each trial 
fit.te <- vector("list", nstudy)
fit.m1 <- vector("list", nstudy)
fit.m2 <- vector("list", nstudy)
fit.y <- vector("list", nstudy)

for (s in 1:nstudy){
  
  # total effects 
  reg.te <- lm(y ~ A + covar1 + covar2, data = studies[[s]])
  
  # TRT -> M1 
  reg.m1 <- lm(M1 ~ A + covar1 + covar2, data = studies[[s]])
  
  # TRT -> M2 
  reg.m2 <-lm(M2 ~ A + covar1 + covar2, data = studies[[s]]) 
  
  # Outcome regression 
  reg.adj.y <- lm(y ~ A*M1*M2 + covar1 + covar2, data = studies[[s]])
  
  fit.te[[s]] <- reg.te
  fit.m1[[s]] <- reg.m1
  fit.m2[[s]] <- reg.m2
  fit.y[[s]] <- reg.adj.y
}

# Part 3: meta-analyze regression analyses --------------------------------

source("run-mvmeta.R")

# set meta-analysis method 
# method = "fixed": fixed-effects meta-analysis
# method = "reml": random-effects meta-analysis

metameth <- "fixed" 

# adjusted outcome regression 
result <- metaBetaCov(fit.y, nstudy, method = metameth)
beta.y <- result[[1]][1,]
cov.y <- result[[2]]

# TRT -> M1  
result <- metaBetaCov(fit.m1, nstudy, method = metameth)
beta.m1 <- result[[1]][1,]
cov.m1 <- result[[2]]

# TRT -> M2
result <- metaBetaCov(fit.m2, nstudy, method = metameth)
beta.m2 <- result[[1]][1,]
cov.m2 <- result[[2]]

# total effects model 
result <- metaBetaCov(fit.te, nstudy, method = metameth)
beta.te <- result[[1]][1,]
cov.te <- result[[2]]

# Part 4: compute decomposed causal effects and variance ------------------

## function to compute causal effects and 
## variance estimates using the delta method 

source("compute-causal-effects.R")

t.meta <- computeCausalEffects(beta.te, beta.m1, beta.m2, beta.y, 
                                 cov.te, cov.m1, cov.m2, cov.y)
t.meta <- data.frame(PointEst = t(t.meta[[1]]), DeltaVar = t(t.meta[[2]]))
colnames(t.meta) <- c("PointEst", "DeltaVar")

## add 95% CIs 
## assuming the effects follow Gaussian distributions
t.meta$lowerCI <- t.meta$PointEst - 1.96*sqrt(t.meta$DeltaVar)
t.meta$upperCI <- t.meta$PointEst + 1.96*sqrt(t.meta$DeltaVar)
 
## t.meta includes the effect estimates and CIs
print(t.meta)


  
