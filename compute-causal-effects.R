## Compute causal mediation effects
## Two mediators, allowing all possible interactions 
## MediatorxMediators
## TreatmentxMediators


computeCausalEffects <- function(beta.te, beta.m1, beta.m2, beta.y, 
                                   cov.te, cov.m1, cov.m2, cov.y){
  cc<-as.matrix(c(0,0))
  
  x10 <- as.numeric(beta.m1[1]+t(cc)%*%beta.m1[-c(1,2)])
  x11 <- as.numeric(beta.m1[1]+beta.m1["A"]+t(cc)%*%beta.m1[-c(1,2)])
  x20 <- as.numeric(beta.m2[1]+t(cc)%*%beta.m2[-c(1,2)])
  x21 <- as.numeric(beta.m2[1]+beta.m2["A"]+t(cc)%*%beta.m2[-c(1,2)])
  
  # Compute the effects of interest
  
  # Total effect: 
  te.model <- beta.te["A"]
  m1 <- 0
  m2 <- 0
  
  # CDE  
  cde <- beta.y["A"] + beta.y["A:M1"]*m1 + beta.y["A:M2"]*m2 + 
    beta.y["A:M1:M2"]*m2*m1
  
  # INT ref M2
  intref_M2 <- (x20-m2)*beta.y["A:M2"]
  
  # INT ref 3-way
  intref_3way <- (x10*x20-m1*m2)*beta.y["A:M1:M2"]
  
  # PNIE
  pnie_M2 <- beta.m2["A"]*beta.y["M2"]
  pnie_M1 <- beta.m1["A"]*beta.y["M1"]
  pnie_3way <- (x11*x21-x10*x20)*beta.y["M1:M2"]
  
  # Med Int 
  medint_M2 <- beta.m2["A"]*beta.y["A:M2"]
  medint_M1 <- beta.m1["A"]*beta.y["A:M1"]
  medint_3way <- (x11*x21-x10*x20)*beta.y["A:M1:M2"]
  
  # Int ref M1 
  intref_M1 <- (x10-m1)*beta.y["A:M1"]
  # alternatively
  # intref_M1 <- te - cde - intref_M2-intref_3way - 
  #   pnie_M2 - pnie_M1 - pnie_3way - medint_M2 - 
  #   medint_M1 - medint_3way
  
  nde.int <- cde + intref_M1 + intref_M2 + intref_3way
  nie.int <- pnie_M2 + pnie_M1 + pnie_3way + medint_M2 + medint_M1 + medint_3way
  
  te <- nde.int + nie.int
  t <- data.frame(cbind(te, nde.int, nie.int, cde, 
                        intref_M2, intref_M1, intref_3way,
                        pnie_M2, pnie_M1, pnie_3way,
                        medint_M2, medint_M1, medint_3way))

  ## delta-method 
  # setting up the point estimate vector and 
  # the variance-covariance diagonal block matrix 
  delta.vec <- c(beta.m1, beta.m2, beta.y) 
  delta.zero <- matrix(0, length(beta.m1), length(beta.m1))
  delta.zero2 <- matrix(0, length(beta.m1), length(beta.y))
  delta.sigma <- rbind(cbind(cov.m1, delta.zero, delta.zero2),
                       cbind(delta.zero, cov.m2, delta.zero2),
                       cbind(t(delta.zero2), t(delta.zero2), cov.y))
  
  
  te.model.var <- cov.te[2,2]
  
  
  # the path-specific causal effects are functions of the regression params
  cde.var <- as.numeric(deltamethod(g = ~ x10, mean = delta.vec, 
                                    cov = delta.sigma, ses = FALSE))
  intref_M1.var <- as.numeric(deltamethod(g = ~ x1*x15, mean = delta.vec, 
                                          cov = delta.sigma, ses = FALSE))
  intref_M2.var <- as.numeric(deltamethod(g = ~ x5*x16, mean = delta.vec, 
                                          cov = delta.sigma, ses = FALSE))
  intref_3way.var <- as.numeric(deltamethod(g = ~ x1*x5*x18, mean = delta.vec, 
                                            cov = delta.sigma, ses = FALSE))
  pnie_M2.var <- as.numeric(deltamethod(g = ~ x6*x12, mean = delta.vec, 
                                        cov = delta.sigma, ses = FALSE))
  pnie_M1.var <- as.numeric(deltamethod(g = ~ x2*x11, mean = delta.vec, 
                                        cov = delta.sigma, ses = FALSE))
  pnie_3way.var <- as.numeric(deltamethod(g = ~ ((x1+x2)*(x5+x6)-(x1*x5))*x17, 
                                          mean = delta.vec, cov = delta.sigma, 
                                          ses = FALSE))
  medint_M2.var <- as.numeric(deltamethod(g = ~ x6*x16, mean = delta.vec, 
                                          cov = delta.sigma, ses = FALSE))
  medint_M1.var <- as.numeric(deltamethod(g = ~ x2*x15, mean = delta.vec, 
                                          cov = delta.sigma, ses = FALSE))
  medint_3way.var <- as.numeric(deltamethod(g = ~ ((x1+x2)*(x5+x6)-(x1*x5))*x18, 
                                            mean = delta.vec, cov = delta.sigma, 
                                            ses = FALSE))
  
  
  # delta method variances for nde.int, nie.int
  # note: te.model here is a placeholder to keep the original order
  
  delta.vec3 <- c(te.model, cde, intref_M2, intref_M1, intref_3way, 
                  pnie_M2, pnie_M1, pnie_3way,
                  medint_M2, medint_M1, medint_3way)
  
  delta.sigma3 <- diag(c(te.model.var, cde.var, intref_M2.var, intref_M1.var, intref_3way.var, 
                         pnie_M2.var, pnie_M1.var, pnie_3way.var,
                         medint_M2.var, medint_M1.var, medint_3way.var))
  
  nde.int.var <- as.numeric(deltamethod(g = ~x2 + x3 + x4 + x5, 
                                        mean = delta.vec3, cov = delta.sigma3, 
                                        ses = FALSE))
  nie.int.var <- as.numeric(deltamethod(g = ~x6 + x7 + x8 + x9 + x10 +x11, 
                                        mean = delta.vec3, cov = delta.sigma3, 
                                        ses = FALSE))
  te.var <- as.numeric(deltamethod(g = ~x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 
                                     x10 + x11, 
                                   mean = delta.vec3, cov = delta.sigma3, 
                                   ses = FALSE))
  
  t.var <- data.frame(cbind(te.var, nde.int.var, nie.int.var,
                            cde.var, intref_M2.var, intref_M1.var, intref_3way.var, 
                            pnie_M2.var, pnie_M1.var, pnie_3way.var,
                            medint_M2.var, medint_M1.var, medint_3way.var))
  
  ## return both causal effect estimates 
  
  return(list(t, t.var))
} 
  
  
  




