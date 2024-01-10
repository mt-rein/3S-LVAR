## this function corrects the standard errors for stepwise estimation
# only for testing!
# load("outputstep1.RData")
# step1step2output <- output1
# load("outputstep3.RData")
# step3output <- output

## auxiliary functions:
# compute rho and epsilon
get_rhoeps <- function(x, param_free, param_est){
  # x:
  #   vector of step1 parameters (lamba, theta, psi)
  #   this is needed for the jacobian() function
  # param_free:
  #   free parameters list from lavInspect
  # param_est:
  #   estimated parameters list from lavInspect
  free_lambdas <- param_free$lambda[param_free$lambda != 0]                   # indices which elements of the lambda matrix are freely estimated
  lambda <- param_est$lambda                                                  # copy the lambda matrix
  lambda[param_free$lambda != 0] <- x[free_lambdas]                           # replace the freely estimated lambdas with the respective elements from x (which is necessary for getting the derivative)
  
  # repeat the above for theta matrix
  free_thetas <- param_free$theta[param_free$theta != 0]
  theta <- param_est$theta
  theta[param_free$theta != 0] <- x[free_thetas]
  
  # repeat the above for psi matrix IF they are not fixed to unity
  free_psis <- param_free$psi[param_free$psi != 0]
  psi <- param_est$psi
  if(length(free_psis) == 0){
    diag(psi) <- 1
  } else {
    psi[param_free$psi != 0] <- x[free_psis]
  }
  
  # obtain sigma and then compute rho and epsilon
  sigma <- lambda %*% psi %*% t(lambda) + theta
  rho <- diag(psi) %*% t(lambda) %*% solve(sigma) %*% lambda
  epsilon <- rho*(1-rho)*diag(psi)
  
  return(c(rho, epsilon))
}

get_LL <- function(x, parameters, data, s2par, s3par){
  n_s2par <- length(s2par)                                                      # number of step2 parameters
  # loop across step2 parameters. Find indices which elements of parameters are 
  # step2 parameters, and replace them with the respective x values (needed to 
  # obtain the derivative)
  for(i in 1:n_s2par){
    row <- which(near(s2par[i], parameters$ustart))
    parameters[row, "est"] <- x[i]
  }
  # do the same for step3 parameters
  for(i in 1:length(s3par)){
    row <- which(parameters$free == i)
    parameters[row, "est"] <- x[i + n_s2par]
  }
  # "fake fit" which is required to obtain the loglikelihood
  fit <- sem(parameters, data = data, missing = "ML", warn = FALSE, do.fit = FALSE)
  LL <- logLik(fit)
  
  return(LL)
}


SEcorrection <- function(step1step2output, step3output){
  # step1step2output:
  #   the object that was generated using the step1step2() function
  # step13output:
  #   the object that was generated using the step3() function
  
  #### 1) Preparations ####
  ## required packages
  library(lavaan)
  library(numDeriv)
  library(dplyr)
  
  ## get information from previous output objects
  fit_step1 <- step1step2output$fit
  fit_step3 <- step3output$fit
  data <- step3output$data
  rho <- step1step2output$rho
  epsilon <- step1step2output$epsilon
  n_step3 <- lavInspect(fit_step3, what = "nobs")
  
  # equation 17 (Bakk et al 2014) gives the corrected sampling variance of the
  # step3 parameters as follows:
  # sigma3corr = sigma3 + solve(H3) %*% C %*% sigma2 %*% t(C) %*% solve(H3)
  # we thus need sigma3, sigma2, H3, and C
  
  #### 2) obtain sigma3 ####
  sigma3 <- vcov(fit_step3)
  
  #### 3) get sigma 2 ####
  # sigma2 is the sampling variance of the step2 parameters (rho and epsilon).
  # this can be obtained from the variances of the step1 parameters (sigma1) and
  # the first-order derivative of the step2 parameters towards the step1
  # parameters.
  
  # obtain the free parameters of step1:
  s1par <- coef(fit_step1, type = "free")
  # numerical derivation using the numDeriv package:
  deriv1 <- jacobian(get_rhoeps, x = s1par,
                       param_free = lavInspect(fit_step1, what = "free"),
                       param_est = lavInspect(fit_step1, what = "est"))
  
  # compute sigma1 and then sigma2:
  sigma1 <- vcov(fit_step1)
  sigma2 <- deriv1 %*% sigma1 %*% t(deriv1)
  
  #### get C ####
  # C is the second-order derivatives of the loglikelihood towards the step3
  # parameters with respect to the step2 parameters. It is a matrix with rows
  # equal to the number of step3 parameters, and columns equal to number of 
  # step2 parameters. 
  # we can obtain it using the hessian() function from numDeriv
  s2par <- c(rho, epsilon)
  s3par <- coef(fit_step3)
  deriv2 <- hessian(get_LL,
                   x = c(s2par, s3par),
                   parameters = parTable(fit_step3),
                   data = data,
                   s2par = s2par,
                   s3par = s3par)
  
  # C is the "lower left" (or upper right) part of the deriv2 matrix
  # (i.e., rows belong to step3 parameters and columns to s2 parameters,
  # or vice versa)
  C <- deriv2[(length(s2par)+1):nrow(deriv2), 1:length(s2par)]
  
  #### get H3 ####
  # H3 contains the second-order derivatives of the loglikelihood towards the
  # step3 parameters.
  # Consequently, it's the lower right part of deriv2 (where rows and columns
  # correspond to the step3 parameters)
  H3 <- deriv2[(length(s2par)+1):nrow(deriv2), (length(s2par)+1):nrow(deriv2)]
  
  ## put everything together and compute the corrected SEs:
  sigma3corr = sigma3 + solve(H3) %*% C %*% sigma2 %*% t(C) %*% solve(H3)
  SEcorr <- sqrt(diag(sigma3corr))
  
  return(SEcorr)
}