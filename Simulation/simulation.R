library(MASS)
library(tidyverse)
library(lavaan)
library(numDeriv)
library(RcppAlgos)


source("sim_VAR.R")
source("step1step2.R")
source("step3.R")
source("SEcorrection.R")
source("auxilliary functions.R")
source("do_sim.R")

# condition grid:
cond <- expand.grid(replication = 1:500,
                    phi_size = c("small", "large"),
                    n = c(25, 50),
                    obs = c(25, 50),
                    rho_gen = c("small", "medium", "large", "very large")
)

# add seeds:
set.seed(123)
cond$seed <- sample(1:nrow(cond)*5, size = nrow(cond), replace = FALSE)
# add iteration number:
cond$iteration <- 1:nrow(cond)                                                  # (unique) number of each iteration


## perform simulation:
output <- lapply(1:nrow(cond), do_sim, cond = cond, outputfile = "output_sim.csv", verbose = TRUE)
# note: the output object is irrelevant, the results are written into the CSV file

# we checked the simulation for errors and warnings, see file "check_warnings.R".
# 516 iterations failed to convergence when using SAM. We re-estimated these 
# using the default bounds when estimating the MM

#### re-estimation with SAM ####
## find iterations where SAM failed to converge
library(tidyverse)
results <- read_csv("Data/output_sim.csv")
failed_iterations <- results$iteration[results$SAM_error]
cond_reestimation <- cond[cond$iteration %in% failed_iterations, ]

## load reestimation function
source("do_sim_reestimation.R")

## run re-estimation
output <- lapply(1:nrow(cond_reestimation), do_sim_reestimation, cond = cond_reestimation, outputfile = "output_sim.csv", verbose = TRUE)