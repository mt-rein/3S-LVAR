# only for testing!
load("examples/outputstep2.RData")
step2output <- output
structuralmodel = NULL

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
  rho <- step2output$rho
  kappa <- step2output$kappa
  fit_step1 <- step2output$fit_step1
  
  ## generate objects for later use
  factors <- lavNames(fit_step1, "lv")                                          # names of latent factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  id <- lavInspect(fit_step1, "cluster")                                        # name of the variable that served as cluster variable in step1
  counts <- data %>% 
    group_by(!!rlang::sym(id)) %>% 
    count()
  max_t <- max(counts$n)
    
  
  #### 2) required data manipulations ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data %>% 
    rename_with(~ factors_ind, all_of(factors))
  
 if(length(factors) == 1){
    data <- data %>% 
      select(all_of(c(id, factors_ind))) %>% 
      group_by(!!rlang::sym(id)) %>% 
      mutate(obs = row_number()) %>% 
      ungroup() %>% 
      pivot_wider(names_from = obs, values_from = all_of(factors_ind), names_prefix = paste0(factors_ind, "_"))
  }
  if(length(factors) > 1){
    data <- data %>% 
      select(all_of(c(id, factors_ind))) %>% 
      group_by(!!rlang::sym(id)) %>% 
      mutate(obs = row_number()) %>% 
      ungroup() %>% 
      pivot_wider(names_from = obs, values_from = all_of(factors_ind))
  }
  
  
  #### 3) write the "measurement" model (factor measured by single indicator) ####
  # (also includes the random intercept)
  measurementmodel <- NULL
  for (fac in factors){
    indicators <- names(data)[grep(fac, names(data))]
    latents <- paste0(fac, "_", 1:length(indicators))
    
    # fix loadings to rho:
    for(col in 1:max_t){
      measurementmodel <- paste0(measurementmodel, latents[col], " =~ ", rho[fac], "*", indicators[col], " \n")
    }
    # fix residual variances to kappa:
    for(col in 1:length(indicators)){
      measurementmodel <- paste0(measurementmodel, indicators[col], " ~~ ", kappa[fac], "*", indicators[col], " \n")
    }
    
    # add population mean (fixed across time):
    measurementmodel <- paste0(measurementmodel, paste0(indicators, collapse = " + "), " ~ grandmean_", fac, "* 1", "\n")
    # add random intercept (i.e., stable deviation of person i from grand mean):
    measurementmodel <- paste0(measurementmodel, paste0("RI_", fac), " =~ ", rho[fac], "*", paste(indicators, collapse = paste0(" + ", rho, "*")), "\n")
    
    measurementmodel <- paste0(measurementmodel, paste0("RI_", fac), " ~~ ", paste0("RI_", fac), "\n")
  }
  
 #### 4) write the structural model ####
  if(is.null(structuralmodel)){
    for (fac1 in factors){
      latents1 <- paste0(fac1, "_", 1:max_t)
      for (fac2 in factors){
        latents2 <- paste0(fac2, "_", 1:max_t)
        # save 
        for(col in 2:max_t){
          # add regression effects:
          structuralmodel <- paste0(structuralmodel, latents1[col], " ~ ", paste0("phi_", fac1, "_", fac2, " * "), latents2[col-1], " \n")
        }
      }
    }
    # add variances at first timepoint:
    
    combs <- comboGeneral(factors, 2, repetition = TRUE)
    for(comb in 1:nrow(combs)){
      for(t in 1:max_t){
        if(t == 1){
          param <- "psi_"
        } else {
          param <- "zeta_"
        }
        structuralmodel <- paste(structuralmodel,
                                 paste0(paste0(combs[comb, 1], "_", t),
                                        paste0(" ~~ ", param, combs[comb, 1], "_", combs[comb, 2], " * "),
                                        paste0(combs[comb, 2], "_", t)),
                                 "\n")
      }
    }
  }
  
  #### 4) combine models and run SEM ####
  fullmodel <- paste(measurementmodel, structuralmodel)
  fit_step3 <- lavaan(fullmodel,
                   data = data,
                   missing = "ML",
                   int.ov.free = FALSE, auto.var = FALSE)
  
  #### 5) extract estimates ####
  ## TO DO: make this more flexible to accomodate for custom SM (provided by user)
  params <- coef(fit_step3)
  phi_names <- unique(names(params))[grep("^phi", unique(names(params)))]
  phis <- purrr::map_dbl(phi_names,\(x) params[names(params) == x][1])
  names(phis) <- phi_names
  
  psi_names <- unique(names(params))[grep("^psi", unique(names(params)))]
  psis <- purrr::map_dbl(psi_names,\(x) params[names(params) == x][1])
  names(psis) <- psi_names
  
  zeta_names <- unique(names(params))[grep("^zeta", unique(names(params)))]
  zetas <- purrr::map_dbl(zeta_names,\(x) params[names(params) == x][1])
  names(zetas) <- zeta_names
  
  estimates <- list("phi" = phis,
                    "psi" = psis,
                    "zeta" = zetas)
  
  
  
  #### 5) build the output ####
  output <- list("fit_step3" = fit_step3,
                 "data" = data,
                 "estimates" = estimates)
  
  return(output)
}