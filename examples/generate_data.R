library(MASS)
library(tidyverse)

source("Simulation/sim_VAR.R")

phi_size <-  "large"
n <- 80
obs <- 50
rho_gen <- "medium"
set.seed(123)

#### 1) set data generation parameters ####
totalvar <- c(1, .3, .3, 1)                                                   # total (co)variance

# AR/CR parameters:
if(phi_size == "small"){
  phi11_pop <- .3                                                             # effect of f1 on f1
  phi12_pop <- .15                                                            # effect of f2 on f1
  phi21_pop <- .15                                                            # effect of f1 on f2
  phi22_pop <- .3                                                             # effect of f2 on f2
  
}

if(phi_size == "large"){
  phi11_pop <- .6                                                             # effect of f1 on f1
  phi12_pop <- .3                                                             # effect of f2 on f1
  phi21_pop <- .3                                                             # effect of f1 on f2
  phi22_pop <- .6                                                             # effect of f2 on f2
}
phimat <- matrix(c(phi11_pop, phi21_pop, phi12_pop, phi22_pop), ncol = 2)     # combine into matrix


# innovation (co)variance matrix:
zetamat <- solve(solve(diag(2*2) - kronecker(phimat, phimat)), totalvar) %>%
  matrix(ncol = 2)                                                            # innovation variance as a function of total variance and phi matrix
zeta1_pop <- zetamat[1, 1]                                                    # innovation variance of x1
zeta2_pop <- zetamat[2, 2]                                                    # innovation variance of x2
zeta12_pop <- zetamat[1, 2]                                                   # innovation covariance between x1 and x2

#### 2) generate factor scores ####
eta <- tibble(id = numeric(),
              obs = numeric(),
              eta1 = numeric(),
              eta2 = numeric())

# loop over all subjects:
for(person in 1:n){
  id <- person                                                                  # ID variable
  mu1 <- rnorm(1, 5, 1)                                                            # person-specific mean on the latent variable
  mu2 <- rnorm(1, 10, 1.5)                                                         # person-specific mean on the latent variable
  
  intercepts <- solve(solve(diag(2)-phimat), c(mu1, mu2))
  
  # generate factor scores for each individual:
  eta_i <- sim_VAR(obs = obs, phi = phimat,
                   zeta = zetamat, intercept = intercepts,
                   burn_in = 0, initintercept = c(mu1, mu2),
                   initvar = matrix(totalvar, ncol = 2))
  
  colnames(eta_i) <- c("eta1", "eta2")
  # add ID variable and observation number to factor scores:
  eta_i <- as_tibble(eta_i) %>%
    mutate(id = !!id,
           obs = row_number())
  
  eta <- full_join(eta, eta_i, by = c("id", "obs", "eta1", "eta2"))
}


#### 3) generate indicators ####
# lambda matrix (loadings), all 1:
loadings <- rep(1, 4)
lambda <- matrix(c(loadings, rep(0, 4), rep(0, 4), loadings),
                 nrow = 8, ncol = 2)
# theta matrix (error variances) depends on rho_gen:
if(rho_gen == "small"){
  errorvar <- rep(3.815, 8)
}
if(rho_gen == "medium"){
  errorvar <- rep(1.603, 8)
}
if(rho_gen == "large"){
  errorvar <- rep(.408, 8)
}
if(rho_gen == "very large"){
  errorvar <- rep(0.0005, 8)
}
theta <- matrix(0, nrow = 8, ncol = 8)                                        # create error covariance matrix
diag(theta) <- errorvar                                                       # impute error variances on theta diagonal (assuming uncorrelated errors)

epsilon <- mvrnorm(nrow(eta), mu = rep(0, 8), Sigma = theta, empirical=T)     # generate errors
# generate indicator scores
# (intercepts left out because we set them to 0):
data <- as.matrix(eta[, c("eta1", "eta2")]) %*% t(lambda) + epsilon %>%
  as.data.frame()
colnames(data) <- c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8")
data$id <- eta$id
data$obs <- eta$obs