# only for testing!
# load("examples/outputstep2.RData")
# step2output <- output
# structuralmodel = NULL

step3 <- function(step2output, structuralmodel = NULL){
  # step2output:
  #   the object that was generated using the step2() function
  # structuralmodel (optional):
  #   a string describing the measurement model using the lavaan syntax.
  #   If ommitted, a VAR model with all auto- and crossregressive parameters
  #   will be specified by default.
  
  #### 1) Preparations ####
  ## required packages:
  library(dplyr)
  library(tidyr)
  library(lavaan)
  library(RcppAlgos)
  
  ## extract objects from step 1 output:
  data <- step2output$data
  lambda_star <- step2output$lambda_star
  theta_star <- step2output$theta_star
  fit_step1 <- step2output$fit_step1
  
  ## generate character vectors for later use
  factors <- lavNames(fit_step1, "lv")                                          # names of latent factors
  factors_lagged <- paste0(factors, "_lag")                                     # names of LAGGED latent factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  factors_ind_lagged <- paste0(factors_ind, "_lag")                             # names of LAGGEd factor score variables (single indicators)
  n_factors <- length(factors)                                                  # numbers of distinct latent constructs
  id <- lavInspect(fit_step1, "cluster")                                        # name of the variable that served as cluster variable in step1
  
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
                            ## fix the loadings to lambda_star (non-lagged and lagged):
                            paste0(factors[fac],
                                   " =~ ",
                                   lambda_star[factors[fac]],
                                   "*",
                                   factors_ind[fac]),
                            "\n",
                            paste0(factors_lagged[fac],
                                   " =~ ",
                                   lambda_star[factors[fac]],
                                   "*",
                                   factors_ind_lagged[fac]),
                            "\n",
                            ## fix the variances to theta_star:
                            paste0(factors_ind[fac],
                                   " ~~ ",
                                   theta_star[factors[fac]],
                                   "*",
                                   factors_ind[fac]),
                            "\n",
                            paste0(factors_ind_lagged[fac],
                                   " ~~ ",
                                   paste(theta_star[factors[fac]]),
                                   "*",
                                   factors_ind_lagged[fac]),
                            "\n"
    )
  }
  
  #### 4) write the structural model ####
  if(is.null(structuralmodel)){
    # get all combinations of DVs (non-lagged, column 1) and IVs (lagged, column 2)
    combinations_phi <- expand.grid(list(factors, factors_lagged))
    # get all combinations of innovation (co)variances
    combinations_zeta <- RcppAlgos::comboGeneral(factors, 2, repetition = TRUE)
    
    ## add phis:
    for(comb in 1:nrow(combinations_phi)){
      structuralmodel <- paste(structuralmodel,
                               paste0(combinations_phi[comb, 1],
                                      " ~ ",
                                      combinations_phi[comb, 2]),
                               "\n"
      )
    }
    
    ## add innovation variances:
    for(comb in 1:nrow(combinations_zeta)){
      structuralmodel <- paste(structuralmodel,
                               paste0(combinations_zeta[comb, 1],
                                      "~~ ",
                                      combinations_zeta[comb, 2]),
                               "\n"
      )
    }
  }
  
  #### 4) combine models and run SEM ####
  fullmodel <- paste(indicatormodel, structuralmodel)
  fit_step3 <- sem(fullmodel,
                   data = data,
                   missing = "ML",
                   cluster = id)
  
  phi <- lavInspect(fit_step3, "est")$beta
  phi <- phi[rowSums(phi) > 0, colSums(phi) > 0]
  
  #### 5) build the output ####
  output <- list("fit_step3" = fit_step3,
                 "data" = data,
                 "phi" = phi)
  
  return(output)
}