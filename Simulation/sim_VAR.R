## this function generates data for a single individual according to a vector autoregressive model
sim_VAR <- function(obs, phi, zeta, intercept, burn_in = 0, initintercept, initvar){
  # obs = number of observations
  # phi = auto-regressive effect (a matrix in case of multiple constructs)
  # zeta = innovation variance (a matrix in case of multiple constructs)
  # intercept = intercept (a vector in case of multiple constructs)
  # burn_in = length of burn in (i.e., data that are generated to remove influence of initial random draw)
  # initintercept = intercept of the first observation (i.e., person mean)
  # initvar = variance of the first observation (i.e., total variance)
  
  
  # create empty dataframe of length obs + burn_in
  data <- as.data.frame(matrix(NA, nrow = burn_in + obs, ncol = 2))
  
  for(i in 1:nrow(data)){
    # simulate the first observation from the person's starting point (= intercept)
    if(i == 1){
      data[i,] <- mvrnorm(n = 1, mu = initintercept, Sigma = initvar)
    }
    
    # then loop through all the rows, predict the current observation from the previous observations, then add random innovation
    if(i > 1){
      predicted <- c(intercept[1] + phi[1,1]*data[i-1, 1] + phi[1,2]*data[i-1, 2], intercept[2] + phi[2,1]*data[i-1, 1] + phi[2,2]*data[i-1, 2])
      data[i, ] <- mvrnorm(n = 1, mu = predicted, Sigma = zeta)
    }
  }
  
  # remove the first rows, depending on length of burn in
  if(burn_in > 0){
    data <- data[-(1:burn_in),] 
  }
  
  rownames(data) <- 1:nrow(data)
  
  return(data)
}