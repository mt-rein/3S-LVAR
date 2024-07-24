correct_nickellsbias <- function(phi, obs){
  k <- ncol(phi) # Number of variables in the VAR model
  I_k <- diag(k) # Identity matrix of order k
  
  # Calculate the bias correction term
  bias <- -(I_k + phi) / (obs - 1)
  corrected_phi <- phi - bias
  return(corrected_phi)
}