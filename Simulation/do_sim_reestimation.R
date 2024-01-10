## define simulation function to re-estimate SAM
# this is identical to the SAM section in do_sim.R, except that default bounds
# in the sam() function are used (instead of mm.args = list(bounds = "none")).
do_sim_reestimation <- function(pos, cond, outputfile, verbose = FALSE){
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished
  replication <- cond$replication[pos]
  iteration <- cond$iteration[pos]
  # get condition levels and set seed:
  phi_size <-  cond$phi_size[pos] %>% as.character()
  n <- cond$n[pos]
  obs <- cond$obs[pos]
  rho_gen <- cond$rho_gen[pos] %>% as.character()
  seed_cond <- cond$seed[pos]
  set.seed(seed_cond)
  
  #### 1) set data generation parameters ####
  totalvar <- c(1, .3, .3, 1)                                                   # total (co)variance
  intercepts <- c(0, 0)                                                         # intercepts have been set to 0
  
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
  
  # start timer for data generation:
  start_datagen <- Sys.time()
  
  #### 2) generate factor scores ####
  eta <- tibble(id = numeric(),
                obs = numeric(),
                eta1 = numeric(),
                eta2 = numeric())
  
  # loop over all subjects:
  for(person in 1:n){
    id <- person                                                                # ID variable
    # generate factor scores for each individual:
    eta_i <- sim_VAR(obs = obs, phi = phimat,
                     zeta = zetamat, intercept = intercepts,
                     burn_in = 0, initintercept = intercepts,
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
  
  # save duration of data generation
  t_datagen <- difftime(Sys.time(), start_datagen, unit = "s")
  
  
  #### 4) SAM ####
  start_SAM <- Sys.time()                                                       # start timer for SAM
  
  data_SEM <- data %>% 
    group_by(id) %>% 
    do(add_row(.)) %>% 
    ungroup() %>% 
    fill(id) %>% 
    group_by(id) %>% 
    mutate(v1lag = dplyr::lag(v1),
           v2lag = dplyr::lag(v2),
           v3lag = dplyr::lag(v3),
           v4lag = dplyr::lag(v4),
           v5lag = dplyr::lag(v5),
           v6lag = dplyr::lag(v6),
           v7lag = dplyr::lag(v7),
           v8lag = dplyr::lag(v8)) %>% 
    ungroup()
  
  
  model_SEM <- "
  f1lag =~ v1lag + v2lag + v3lag + v4lag
  f2lag =~ v5lag + v6lag + v7lag + v8lag
  f1 =~ v1 + v2 + v3 + v4
  f2 =~ v5 + v6 + v7 + v8
  
  f1 ~ phi_f1_f1_lag*f1lag + phi_f1_f2_lag*f2lag
  f2 ~ phi_f2_f1_lag*f1lag + phi_f2_f2_lag*f2lag
  f1 ~~ zeta_f1_f1*f1
  f2 ~~ zeta_f2_f2*f2
  f1 ~~ zeta_f1_f2*f2
  "
  
  SAM <- run_SAM(model_SEM, data = data_SEM, missing = "ML",
                 sam.method = "local")
  # extract error/warning messages (if applicable):
  SAM_warning <- ifelse(is_empty(SAM$warnings),
                        FALSE, TRUE)
  SAM_warning_text <- ifelse(is_empty(SAM$warnings),
                             "",
                             paste(c(SAM$warnings),
                                   collapse = "; ")
  )
  SAM_error <- ifelse(is_empty(SAM$result$error),
                      FALSE, TRUE)
  SAM_error_text <- ifelse(is_empty(SAM$result$error),
                           "",
                           paste(c(SAM$result$error),
                                 collapse = "; ")
  )
  
  # if SAM was succesful, extract results:
  if(is_empty(SAM$result$error)){
    # save duration of SAM:
    output_SAM <- SAM$result$result
    t_SAM <- difftime(Sys.time(), start_SAM, unit = "s")                        # save duration for SAM
    
    ## extract results:
    # regression parameters and innovation variances
    SAM_parameters <- parTable(output_SAM)
    SAM_phi11 <- SAM_parameters %>% 
      filter(label == "phi_f1_f1_lag") %>% 
      dplyr::select(est) %>% as.numeric()
    SAM_phi12 <- SAM_parameters %>% 
      filter(label == "phi_f1_f2_lag") %>% 
      dplyr::select(est) %>% as.numeric()
    SAM_phi21 <- SAM_parameters %>% 
      filter(label == "phi_f2_f1_lag") %>% 
      dplyr::select(est) %>% as.numeric()
    SAM_phi22 <- SAM_parameters %>% 
      filter(label == "phi_f2_f2_lag") %>% 
      dplyr::select(est) %>% as.numeric()
    
    SAM_zeta1 <- SAM_parameters %>%
      filter(label == "zeta_f1_f1") %>% 
      dplyr::select(est) %>% as.numeric()
    SAM_zeta2 <- SAM_parameters %>%
      filter(label == "zeta_f2_f2") %>% 
      dplyr::select(est) %>% as.numeric()
    SAM_zeta12 <- SAM_parameters %>%
      filter(label == "zeta_f1_f2") %>% 
      dplyr::select(est) %>% as.numeric()
    
    # standard errors
    SAM_phi11_se <- SAM_parameters %>% 
      filter(label == "phi_f1_f1_lag") %>% 
      dplyr::select(se) %>% as.numeric()
    SAM_phi12_se <- SAM_parameters %>% 
      filter(label == "phi_f1_f2_lag") %>% 
      dplyr::select(se) %>% as.numeric()
    SAM_phi21_se <- SAM_parameters %>% 
      filter(label == "phi_f2_f1_lag") %>% 
      dplyr::select(se) %>% as.numeric()
    SAM_phi22_se <- SAM_parameters %>% 
      filter(label == "phi_f2_f2_lag") %>% 
      dplyr::select(se) %>% as.numeric()
    
    SAM_zeta1_se <- SAM_parameters %>%
      filter(label == "zeta_f1_f1") %>% 
      dplyr::select(se) %>% as.numeric()
    SAM_zeta2_se <- SAM_parameters %>%
      filter(label == "zeta_f2_f2") %>% 
      dplyr::select(se) %>% as.numeric()
    SAM_zeta12_se <- SAM_parameters %>%
      filter(label == "zeta_f1_f2") %>% 
      dplyr::select(se) %>% as.numeric()
  } else{
    # if SAM was not succesful, set results to NA
    t_SAM <- NA
    
    SAM_phi11 <- NA
    SAM_phi12 <- NA
    SAM_phi21 <- NA
    SAM_phi22 <- NA
    
    SAM_zeta1 <- NA
    SAM_zeta2 <- NA
    SAM_zeta12 <- NA
    
    SAM_phi11_se <- NA
    SAM_phi12_se <- NA
    SAM_phi21_se <- NA
    SAM_phi22_se <- NA
    
    SAM_zeta1_se <- NA
    SAM_zeta2_se <- NA
    SAM_zeta12_se <- NA
  }
  
  
  # compile output:
  output <- c("iteration" = iteration, "replication" = replication,
              "phi_size" = phi_size, "n" = n, "obs" = obs, "rho_gen" = rho_gen,
              "t_datagen" = t_datagen, "t_SAM" = t_SAM,
              "phi11_pop" = phi11_pop, "phi12_pop" = phi12_pop, "phi21_pop" = phi21_pop, "phi22_pop" = phi22_pop,
              "zeta1_pop" = zeta1_pop, "zeta2_pop" = zeta2_pop, "zeta12_pop" = zeta12_pop,
              "SAM_phi11" = SAM_phi11, "SAM_phi12" = SAM_phi12,
              "SAM_phi21" = SAM_phi21, "SAM_phi22" = SAM_phi22,
              "SAM_zeta1" = SAM_zeta1, "SAM_zeta2" = SAM_zeta2, "SAM_zeta12" = SAM_zeta12,
              "SAM_phi11_se" = SAM_phi11_se, "SAM_phi12_se" = SAM_phi12_se,
              "SAM_phi21_se" = SAM_phi21_se, "SAM_phi22_se" = SAM_phi22_se,
              "SAM_zeta1_se" = SAM_zeta1_se, "SAM_zeta2_se" = SAM_zeta2_se, "SAM_zeta12_se" = SAM_zeta12_se,
              "SAM_warning" = SAM_warning,"SAM_error" = SAM_error,
              "seed" = seed_cond, "pos" = pos,
              "SAM_warning_text" = SAM_warning_text, "SAM_error_text" = SAM_error_text)
  
  for(i in 34:35){
    output[i] <- str_squish(output[i])                                          # removes all whitespace and linebreaks from the error and warning strings
    output[i] <- gsub(",", "", output[i])                                       # removes all commata from error and warning strings (to prevent messing up the CSV file)
  }
  
  
  if(pos == 1){
    write.table(t(output), file = outputfile, append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  }
  
  if(pos > 1){
    write.table(t(output), file = outputfile, append = TRUE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  }
  
  if(verbose == TRUE){
    print(paste("Simulation", pos, "completed at", Sys.time())) # prints when a replication is done (as a sign that R did not crash)
  }
  
  invisible(output)
}