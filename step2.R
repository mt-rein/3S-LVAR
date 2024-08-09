# only for testing!
# load("examples/outputstep1.RData")
# step1output <- output

step2 <- function(step1output){
  # step1output:
  #   the object that was generated using the step1() function
  
  ## Preparations
  library(lavaan)
  fit_step1 <- step1output$fit_step1
  data <- step1output$data
  
  ## compute factor scores:
  factorscores <- lavPredict(fit_step1, assemble = TRUE)
  # append factor scores to original data:
  data <- cbind(data, factorscores)
  
  ## compute lambda_star and theta_star
  # These values will be used in step 3 to account for measurement error
  psi <- lavInspect(fit_step1, what = "est")$psi
  lambda <- lavInspect(fit_step1, what = "est")$lambda
  sigma <- fitted(fit_step1)$cov
  lambda_star <- diag(psi %*% t(lambda) %*% solve(sigma) %*% lambda)
  theta_star <- lambda_star*(1-lambda_star)*diag(psi)
  
  # assemble output
  output <- list("data" = data,
                 "lambda_star" = lambda_star,
                 "theta_star" = theta_star,
                 "fit_step1" = fit_step1)
  return(output)
}