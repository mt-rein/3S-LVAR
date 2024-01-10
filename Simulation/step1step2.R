## this function performs step1 and step2 of 3S-LVAR
# only for testing!
# data <- read_csv("exampledata.csv")
# measurementmodel <- "
#   f1 =~ V1 + V2 + V3 + V4 + V5
#   f2 =~ V6 + V7 + V8 + V9 + V10
# 
#   V1 ~ 0*1
#   V2 ~ 0*1
#   V3 ~ 0*1
#   V4 ~ 0*1
#   V5 ~ 0*1
#   V6 ~ 0*1
#   V7 ~ 0*1
#   V8 ~ 0*1
#   V9 ~ 0*1
#   V10 ~ 0*1
#   "
# id = "id"

step1step2 <- function(data, measurementmodel, id){
  # data:
  #   a data frame with the indicator variables, as well as id and group (for
  #   step 2) variables
  # measurementmodel:
  #   a string describing the measurement model using the lavaan syntax
  # id:
  #   a character that indicates the id variable (the variable that indicates
  #   which observations belong to which person)
  
  ## Preparations
  library(lavaan)
  indicators <- lavNames(measurementmodel, "ov.ind")                            # vector of indicator variable names
  factors <- lavNames(measurementmodel, "lv")                                   # vector of latent factor names
  
  # estimate the measurement model:
  fit <- cfa(measurementmodel,
             data = data,
             orthogonal = TRUE,
             missing = "ML",
             cluster = id)

  # compute factor scores:
  factorscores <- lavPredict(fit, assemble = TRUE, append.data = TRUE)
  # append factor scores to original data:
  data <- cbind(data, factorscores[, factors])
  
  # compute rho:
  # (rho indicates the model-based reliability.
  # These values will be used in step 2 to account for measurement error.)
  psi <- lavInspect(fit, what = "est")$psi
  lambda <- lavInspect(fit, what = "est")$lambda
  sigma <- fitted(fit)$cov
  rho <- diag(psi %*% t(lambda) %*% solve(sigma) %*% lambda)
  
  # compute epsilon:
  # (epsilon indicates the uncertainty of the factor score estimation)
  epsilon <- rho*(1-rho)*diag(psi)
  
  # assemble output
  output <- list("data" = data,
                 "rho" = rho,
                 "epsilon" = epsilon,
                 "fit" = fit)
  return(output)
}