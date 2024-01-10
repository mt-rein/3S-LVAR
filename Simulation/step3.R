## This function performs step 3 of 3S-LVAR
# only for testing!
# load("outputstep1.RData")
# step1step2output <- output1
# id = "id"

step3 <- function(step1step2output, id){
  # step1step2output:
  #   the object that was generated using the step1step2() function
  # id:
  #   a character that indicates the id variable (the variable that indicates
  #   which observations belong to which person)
  # group:
  #   a character that indicates the grouping variable (the variable that
  #   indicates which person belongs to which higher level group)
  
  #### 1) Preparations ####
  ## required packages:
  library(dplyr)
  library(tidyr)
  library(lavaan)
  library(RcppAlgos)
  
  ## extract objects from step 1 output:
  data <- step1step2output$data
  rho <- step1step2output$rho
  epsilon <- step1step2output$epsilon
  fit_step1 <- step1step2output$fit
  
  ## generate character vectors for later use
  factors <- names(rho)                                                         # names of latent factors
  factors_lagged <- paste0(factors, "_lag")                                     # names of LAGGED latent factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  factors_ind_lagged <- paste0(factors_ind, "_lag")                             # names of LAGGEd factor score variables (single indicators)
  n_factors <- length(factors)                                                  # numbers of distinct latent constructs
  
  #### 2) required data manipulations ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data %>% 
    rename_with(~ factors_ind, all_of(factors))
  
  ## add additional row per individual and generate lagged variables
  data <- data %>% 
    group_by(!!rlang::sym(id)) %>%
    do(add_row(.)) %>% 
    ungroup() %>% 
    fill(!!rlang::sym(id))                                                      # this fills in the id variable for the new rows
  
  for(var in factors_ind){
    data <- data %>% 
      group_by(!!rlang::sym(id)) %>% 
      mutate("{var}_lag" := dplyr::lag(!!rlang::sym(var))) %>%                  # create the lagged variables with dynamic variable names
      ungroup()
  }
  
  #### 3) write the "measurement" model (factor measured by single indicator) ####
  indicatormodel <- NULL
  for (fac in 1:n_factors){
    indicatormodel <- paste(indicatormodel,
                            ## fix the loadings to rho (non-lagged and lagged):
                            paste0(factors[fac],
                                   " =~ ",
                                   rho[factors[fac]],
                                   "*",
                                   factors_ind[fac]),
                            "\n",
                            paste0(factors_lagged[fac],
                                   " =~ ",
                                   rho[factors[fac]],
                                   "*",
                                   factors_ind_lagged[fac]),
                            "\n",
                            ## fix the variances to epsilon:
                            paste0(factors_ind[fac],
                                   " ~~ ",
                                   epsilon[factors[fac]],
                                   "*",
                                   factors_ind[fac]),
                            "\n",
                            paste0(factors_ind_lagged[fac],
                                   " ~~ ",
                                   paste(epsilon[factors[fac]]),
                                   "*",
                                   factors_ind_lagged[fac]),
                            "\n"
    )
  }
  
  #### 4) write the structural model ####
  # get all combinations of DVs (non-lagged, column 1) and IVs (lagged, column 2)
  combinations_phi <- expand.grid(list(factors, factors_lagged))
  # get all combinations of innovation (co)variances
  combinations_zeta <- RcppAlgos::comboGeneral(factors, 2, repetition = TRUE)
  
  structuralmodel <- NULL
  ## add phis:
  for(comb in 1:nrow(combinations_phi)){
    structuralmodel <- paste(structuralmodel,
                             paste0(combinations_phi[comb, 1],
                                    " ~ ",
                                    paste0("phi_",
                                           combinations_phi[comb, 1], "_",
                                           combinations_phi[comb, 2]),
                                    "*",
                                    combinations_phi[comb, 2]),
                             "\n"
    )
  }
  # Note: for example, "phi_A_B_lag" means that this parameter refers to the
  # regression of A on B_lag
  
  ## add innovation variances:
  for(comb in 1:nrow(combinations_zeta)){
    structuralmodel <- paste(structuralmodel,
                             paste0(combinations_zeta[comb, 1],
                                    "~~ ",
                                    paste0("zeta_",
                                           combinations_zeta[comb, 1], "_",
                                           combinations_zeta[comb, 2]),
                                    "*",
                                    combinations_zeta[comb, 2]),
                             "\n"
    )
  }
  
  #### 4) combine models and run SEM ####
  fullmodel <- paste(indicatormodel, structuralmodel)
  fit_step3 <- sem(fullmodel,
                   data = data,
                   missing = "ML",
                   cluster = id)
  
  
  #### 5) build the output ####
  output <- list("beta" = lavInspect(fit_step3, "est")$beta,
                 "all_estimates" = lavInspect(fit_step3, "est"),
                 "standarderrors" = lavInspect(fit_step3, "se"),
                 "fit" = fit_step3,
                 "data" = data)
  
  return(output)
}