## which (datagenerating) values for theta lead to a (datagenerating) value for rho that is in line with the levels of the simulation?
lambda <- matrix(c(rep(1, 4), rep(0, 4), rep(0, 4), rep(1, 4)), nrow = 8, ncol = 2)
psi <- matrix(c(1, .3, .3, 1), nrow = 2)
theta <- matrix(0, nrow = 8, ncol = 8) # create error covariance matrix
diag(theta) <- 0.0005
sigma <- lambda %*% psi %*% t(lambda) + theta

rho <- diag(psi %*% t(lambda) %*% solve(sigma) %*% lambda)
epsilon <- rho*(1-rho)*diag(psi)


# theta: 3.815 --> rho .5
# theta: 1.603 --> rho .7
# theta: .408 --> rho .9
# theta: 0.035 --> rho .99
# theta: 0.0005 --> rho .999