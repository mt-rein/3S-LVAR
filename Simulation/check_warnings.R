## this script checks the initial simulation results on errors and warning messages

library(tidyverse)

results <- read_csv("Data/output_sim.csv") %>% 
  mutate(phi_size = factor(phi_size, levels = c("small", "large")),
         rho_gen = factor(rho_gen, levels = c("small", "medium", "large", "very large"))
  )
# load data and turn two variables into factors

#### errors ####
sum(results$LVAR_step1step2_error)
sum(results$LVAR_step3_error)
sum(results$SEcorr_error)
sum(results$NFS_error)
sum(results$SAM_error)
sum(results$SEM_error)
# 516 errors for SAM, 0 for all others

unique(results$SAM_error_text)
# "number of items to replace is not a multiple of replacement length; Sigma.11[par.idx par.idx] <- sigma.11[keep.idx keep.idx drop = FALSE]"
unique(results$SAM_warning_text[results$SAM_error])
# these errors seem to be caused by non-convergence of the measurement model

#### missing results ####
results %>% 
  summarise(across(21:83, ~ sum(is.na(.x)))) %>% 
  print(width = Inf)
# 516 for SAM --> due do the errors/non-convergence

#### warnings ####
# LVAR step1/step2:
sum(results$LVAR_step1step2_warning)
# 1002
head(unique(results$LVAR_step1step2_warning_text))
sum(grepl("positive definite", results$LVAR_step1step2_warning_text,  fixed = TRUE))
# all of these warnings are about the variance-covariance matrix of the estimated 
# parameters not being positive definite. We found that these warnings are caused 
# by using the cluster-robust standard  errors. The reason is that the "sandwich" 
# matrix is not  always full rank. We chose to ignore the warning, however, after 
# confirming with other software (LatentGold) that the results can be trusted.

# LVAR step3:
sum(results$LVAR_step3_warning)
# 79
head(unique(results$LVAR_step3_warning_text))
sum(grepl("positive definite", results$LVAR_step3_warning_text,  fixed = TRUE))
# same as for step1/step2 ("vcov not positive definite" warning)

# SE correction:
sum(results$SEcorr_warning)
# 0
# no warnings for the SE correction

# NFS:
sum(results$NFS_warning)
# 16000
unique(results$NFS_warning_text)
# For NFS, lavaan warned about the deletion of missing variables. This is caused 
# by creating the lagged variables and not problematic

# SAM:
sum(results$SAM_warning)
# 516
unique(results$SAM_warning_text)
# all warnings due to non-convergence (see above)

# SEM:
sum(results$SEM_warning)
# 16000
head(unique(results$SEM_warning_text))
sum(grepl("positive definite", results$SEM_warning_text,  fixed = TRUE))
# same as for LVAR ("vcov not positive definite" warning)



#### SAM re-estimation ####
results_reestimation <- read_csv("Data/output_sim_reestimation.csv") %>% 
  mutate(phi_size = factor(phi_size, levels = c("small", "large")),
         rho_gen = factor(rho_gen, levels = c("small", "medium", "large", "very large"))
  )

sum(results_reestimation$SAM_error)
# 153 errors (out of 516)
unique(results$SAM_warning_text)
sum(results_reestimation$SAM_warning)
# all warnings due to non-convergence (see above)
# --> convergence was improved, but 153 still failed to converge