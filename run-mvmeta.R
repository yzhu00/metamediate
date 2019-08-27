

#### Meta-analysis ####

# metaBetaCov
# a function that gets meta-analyzed coefs and vcovs for a given model
# mdl.list: a list of regression models from individual studies 
# nstudy: number of studies 
# method = "fixed": fixed-effects meta-analysis
# method = "reml": random-effects meta-analysis

metaBetaCov <- function(mdl.list, nstudy, method){
  nbeta <- length(mdl.list[[1]]$coeff)
  
  # save beta coef and vcov from each outcome model 
  beta_study <- matrix(vector(), nbeta, 0)
  Sigma <- array(0, c(nbeta, nbeta, nstudy))
  vec <- matrix(0,nbeta*(nbeta+1)/2,nstudy)
  
  for (i in 1:nstudy){
    beta_study <- cbind(beta_study, mdl.list[[i]]$coeff)
    Sigma[, , i] <-  vcov(mdl.list[[i]]) 
    vec[,i] <- vechMat(Sigma[,,i], diag = TRUE)
  }
  
  o <- t(beta_study)
  X <- matrix(1,nstudy,1)
  S <- t(vec)
  
  Id = diag(nstudy)
  
  Sigma_big<-bdiag(lapply(1:nstudy, function(x){Sigma[,,x]}))
  W = diag(1,nbeta,nbeta)
  for (w in 2:nstudy){
    W = rbind(W, diag(1,nbeta,nbeta))
  }
  
  
  ### Method = fixed 
  if (method == "fixed"){
    model_fixed <- mvmeta.fit(X,o,S, method="fixed")
    meta.beta <- model_fixed$coef
    meta.cov <- model_fixed$vcov
  }
  
  ### Method = reml 
  if (method == "reml"){
    model_reml <- mvmeta.fit(X,o,S, method="reml",control=list(showiter=TRUE,igls.iter=10, maxiter = 500)) 
    meta.beta <- model_reml$coef
    
    Tau_big<-kronecker(Id,model_reml$Psi)
    Omega<-(Sigma_big+Tau_big)
    Omegainv = solve(Omega)                  
    meta.cov <- as.matrix(solve(t(W)%*%Omegainv%*%W))
  }

  
  return(list(beta = meta.beta, cov = meta.cov))
}

