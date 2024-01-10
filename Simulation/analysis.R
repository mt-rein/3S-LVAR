library(tidyverse)
# directly load the cleaned data instead of creating them:
# results_clean <- read_csv("Data/results_clean.csv")


## load data and merge original and re-estimation data set
results_original <- read_csv("Data/output_sim.csv") %>% 
  mutate(phi_size = factor(phi_size, levels = c("small", "large")),
         rho_gen = factor(rho_gen, levels = c("small", "medium", "large", "very large"))
  ) %>% 
  arrange(iteration)                                                            # sort by iteration
results_reestimation <- read_csv("Data/output_sim_reestimation.csv") %>% 
  mutate(phi_size = factor(phi_size, levels = c("small", "large")),
         rho_gen = factor(rho_gen, levels = c("small", "medium", "large", "very large"))
  ) %>% 
  arrange(iteration)                                                            # sort by iteration
# load data and turn two variables into factors

# update the "original" results where SAM failed with the re-estimation results
results <- results_original
results[results$SAM_error, names(results_reestimation)] <- results_reestimation

# remove 153 iterations where SAM failed to converge in both settings (with both default and no bounds)
# (see check_warnings.R file)
results_clean <- results %>% 
  filter(!SAM_error)

# check if any parameters were not estimated
results_clean %>% 
  summarise(across(21:83, ~ sum(is.na(.)))) %>% 
  print(width = Inf)

# how many replications remain in each condition?
results_clean %>% 
  group_by(phi_size, n, obs, rho_gen) %>% 
  summarise(n_rep = n()) %>% 
  arrange(n_rep)
# largest number of removed replications for a condition is 73


## time
# overall
results_clean %>% 
  summarise(m_time_datagen = mean(t_datagen, na.rm = TRUE),
            sd_time_datagen = sd(t_datagen, na.rm = TRUE),
            m_time_step1step2 = mean(t_step1step2, na.rm = TRUE),
            sd_time_step1step2 = sd(t_step1step2, na.rm = TRUE),
            m_time_step3 = mean(t_step3, na.rm = TRUE),
            sd_time_step3 = sd(t_step3, na.rm = TRUE),
            m_time_SEcorr = mean(t_SEcorr, na.rm = TRUE),
            sd_time_SEcorr = sd(t_SEcorr, na.rm = TRUE),
            m_time_NFS = mean(t_NFS, na.rm = TRUE),
            sd_time_NFS = sd(t_NFS, na.rm = TRUE),
            m_time_SAM = mean(t_SAM, na.rm = TRUE),
            sd_time_SAM = sd(t_SAM, na.rm = TRUE),
            m_time_SEM = mean(t_SEM, na.rm = TRUE),
            sd_time_SEM = sd(t_SEM, na.rm = TRUE)
            ) %>% 
  print(width = Inf)

# by phi size
results_clean %>% 
  group_by(phi_size) %>% 
  summarise(m_time_datagen = mean(t_datagen, na.rm = TRUE),
            sd_time_datagen = sd(t_datagen, na.rm = TRUE),
            m_time_step1step2 = mean(t_step1step2, na.rm = TRUE),
            sd_time_step1step2 = sd(t_step1step2, na.rm = TRUE),
            m_time_step3 = mean(t_step3, na.rm = TRUE),
            sd_time_step3 = sd(t_step3, na.rm = TRUE),
            m_time_SEcorr = mean(t_SEcorr, na.rm = TRUE),
            sd_time_SEcorr = sd(t_SEcorr, na.rm = TRUE),
            m_time_NFS = mean(t_NFS, na.rm = TRUE),
            sd_time_NFS = sd(t_NFS, na.rm = TRUE),
            m_time_SAM = mean(t_SAM, na.rm = TRUE),
            sd_time_SAM = sd(t_SAM, na.rm = TRUE),
            m_time_SEM = mean(t_SEM, na.rm = TRUE),
            sd_time_SEM = sd(t_SEM, na.rm = TRUE)
  ) %>% 
  print(width = Inf)
# no differences

## by n
results_clean %>% 
  group_by(n) %>% 
  summarise(m_time_datagen = mean(t_datagen, na.rm = TRUE),
            sd_time_datagen = sd(t_datagen, na.rm = TRUE),
            m_time_step1step2 = mean(t_step1step2, na.rm = TRUE),
            sd_time_step1step2 = sd(t_step1step2, na.rm = TRUE),
            m_time_step3 = mean(t_step3, na.rm = TRUE),
            sd_time_step3 = sd(t_step3, na.rm = TRUE),
            m_time_SEcorr = mean(t_SEcorr, na.rm = TRUE),
            sd_time_SEcorr = sd(t_SEcorr, na.rm = TRUE),
            m_time_NFS = mean(t_NFS, na.rm = TRUE),
            sd_time_NFS = sd(t_NFS, na.rm = TRUE),
            m_time_SAM = mean(t_SAM, na.rm = TRUE),
            sd_time_SAM = sd(t_SAM, na.rm = TRUE),
            m_time_SEM = mean(t_SEM, na.rm = TRUE),
            sd_time_SEM = sd(t_SEM, na.rm = TRUE)
  ) %>% 
  print(width = Inf)
# as expected, data generation takes twice as long
# no (real) differences otherwise

## by obs
results_clean %>% 
  group_by(obs) %>% 
  summarise(m_time_datagen = mean(t_datagen, na.rm = TRUE),
            sd_time_datagen = sd(t_datagen, na.rm = TRUE),
            m_time_step1step2 = mean(t_step1step2, na.rm = TRUE),
            sd_time_step1step2 = sd(t_step1step2, na.rm = TRUE),
            m_time_step3 = mean(t_step3, na.rm = TRUE),
            sd_time_step3 = sd(t_step3, na.rm = TRUE),
            m_time_SEcorr = mean(t_SEcorr, na.rm = TRUE),
            sd_time_SEcorr = sd(t_SEcorr, na.rm = TRUE),
            m_time_NFS = mean(t_NFS, na.rm = TRUE),
            sd_time_NFS = sd(t_NFS, na.rm = TRUE),
            m_time_SAM = mean(t_SAM, na.rm = TRUE),
            sd_time_SAM = sd(t_SAM, na.rm = TRUE),
            m_time_SEM = mean(t_SEM, na.rm = TRUE),
            sd_time_SEM = sd(t_SEM, na.rm = TRUE)
  ) %>% 
  print(width = Inf)
# same as when n increases (as to be expected)

## by rho
results_clean %>% 
  group_by(rho_gen) %>% 
  summarise(m_time_datagen = mean(t_datagen, na.rm = TRUE),
            sd_time_datagen = sd(t_datagen, na.rm = TRUE),
            m_time_step1step2 = mean(t_step1step2, na.rm = TRUE),
            sd_time_step1step2 = sd(t_step1step2, na.rm = TRUE),
            m_time_step3 = mean(t_step3, na.rm = TRUE),
            sd_time_step3 = sd(t_step3, na.rm = TRUE),
            m_time_SEcorr = mean(t_SEcorr, na.rm = TRUE),
            sd_time_SEcorr = sd(t_SEcorr, na.rm = TRUE),
            m_time_NFS = mean(t_NFS, na.rm = TRUE),
            sd_time_NFS = sd(t_NFS, na.rm = TRUE),
            m_time_SAM = mean(t_SAM, na.rm = TRUE),
            sd_time_SAM = sd(t_SAM, na.rm = TRUE),
            m_time_SEM = mean(t_SEM, na.rm = TRUE),
            sd_time_SEM = sd(t_SEM, na.rm = TRUE)
  ) %>% 
  print(width = Inf)
# no differences except for very large rho (= .999)
# in this condition, 3S-LVAR step 1/2, SEM and SAM take about twice as long

#### performance in each condition ####
performance <- results_clean %>%
  group_by(rho_gen, obs, n, phi_size) %>% 
  summarise(
    ## bias:
    AB_phi11_LVAR = mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_LVAR = (mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_LVAR = mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_LVAR = (mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_LVAR = mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_LVAR = (mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_LVAR = mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_LVAR = (mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_LVAR = mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_LVAR = (mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_LVAR = mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_LVAR = (mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_LVAR = mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_LVAR = (mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # LVAR
    AB_phi11_NFS = mean(NFS_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_NFS = (mean(NFS_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_NFS = mean(NFS_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_NFS = (mean(NFS_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_NFS = mean(NFS_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_NFS = (mean(NFS_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_NFS = mean(NFS_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_NFS = (mean(NFS_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_NFS = mean(NFS_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_NFS = (mean(NFS_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_NFS = mean(NFS_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_NFS = (mean(NFS_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_NFS = mean(NFS_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_NFS = (mean(NFS_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # NFS
    AB_phi11_SAM = mean(SAM_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SAM = (mean(SAM_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SAM = mean(SAM_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SAM = (mean(SAM_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_SAM = mean(SAM_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_SAM = (mean(SAM_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_SAM = mean(SAM_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_SAM = (mean(SAM_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_SAM = mean(SAM_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_SAM = (mean(SAM_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_SAM = mean(SAM_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_SAM = (mean(SAM_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_SAM = mean(SAM_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_SAM = (mean(SAM_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # SAM
    AB_phi11_SEM = mean(SEM_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SEM = (mean(SEM_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SEM = mean(SEM_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SEM = (mean(SEM_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_SEM = mean(SEM_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_SEM = (mean(SEM_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_SEM = mean(SEM_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_SEM = (mean(SEM_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_SEM = mean(SEM_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_SEM = (mean(SEM_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_SEM = mean(SEM_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_SEM = (mean(SEM_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_SEM = mean(SEM_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_SEM = (mean(SEM_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # SEM
    ## RMSE:
    RMSE_phi11_LVAR = mean(sqrt((LVAR_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_LVAR = mean(sqrt((LVAR_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_LVAR = mean(sqrt((LVAR_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_LVAR = mean(sqrt((LVAR_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_LVAR = mean(sqrt((LVAR_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_LVAR = mean(sqrt((LVAR_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_LVAR = mean(sqrt((LVAR_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # LVAR
    RMSE_phi11_NFS = mean(sqrt((NFS_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_NFS = mean(sqrt((NFS_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_NFS = mean(sqrt((NFS_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_NFS = mean(sqrt((NFS_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_NFS = mean(sqrt((NFS_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_NFS = mean(sqrt((NFS_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_NFS = mean(sqrt((NFS_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # NFS
    RMSE_phi11_SAM = mean(sqrt((SAM_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SAM = mean(sqrt((SAM_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_SAM = mean(sqrt((SAM_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_SAM = mean(sqrt((SAM_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_SAM = mean(sqrt((SAM_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_SAM = mean(sqrt((SAM_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_SAM = mean(sqrt((SAM_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # SAM
    RMSE_phi11_SEM = mean(sqrt((SEM_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SEM = mean(sqrt((SEM_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_SEM = mean(sqrt((SEM_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_SEM = mean(sqrt((SEM_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_SEM = mean(sqrt((SEM_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_SEM = mean(sqrt((SEM_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_SEM = mean(sqrt((SEM_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # SEM
    ## Standard error recover:
    SER_phi11_LVAR = mean(LVAR_phi11_se, na.rm = TRUE) / sd(LVAR_phi11, na.rm = TRUE),
    SER_phi22_LVAR = mean(LVAR_phi22_se, na.rm = TRUE) / sd(LVAR_phi22, na.rm = TRUE),
    SER_phi12_LVAR = mean(LVAR_phi12_se, na.rm = TRUE) / sd(LVAR_phi12, na.rm = TRUE),
    SER_phi21_LVAR = mean(LVAR_phi21_se, na.rm = TRUE) / sd(LVAR_phi21, na.rm = TRUE),
    SER_zeta1_LVAR = mean(LVAR_zeta1_se, na.rm = TRUE) / sd(LVAR_zeta1, na.rm = TRUE),
    SER_zeta2_LVAR = mean(LVAR_zeta2_se, na.rm = TRUE) / sd(LVAR_zeta2, na.rm = TRUE),
    SER_zeta12_LVAR = mean(LVAR_zeta12_se, na.rm = TRUE) / sd(LVAR_zeta12, na.rm = TRUE),
    SER_phi11corr_LVAR = mean(LVAR_phi11_secorr, na.rm = TRUE) / sd(LVAR_phi11, na.rm = TRUE),
    SER_phi22corr_LVAR = mean(LVAR_phi22_secorr, na.rm = TRUE) / sd(LVAR_phi22, na.rm = TRUE),
    SER_phi12corr_LVAR = mean(LVAR_phi12_secorr, na.rm = TRUE) / sd(LVAR_phi12, na.rm = TRUE),
    SER_phi21corr_LVAR = mean(LVAR_phi21_secorr, na.rm = TRUE) / sd(LVAR_phi21, na.rm = TRUE),
    SER_zeta1corr_LVAR = mean(LVAR_zeta1_secorr, na.rm = TRUE) / sd(LVAR_zeta1, na.rm = TRUE),
    SER_zeta2corr_LVAR = mean(LVAR_zeta2_secorr, na.rm = TRUE) / sd(LVAR_zeta2, na.rm = TRUE),
    SER_zeta12corr_LVAR = mean(LVAR_zeta12_secorr, na.rm = TRUE) / sd(LVAR_zeta12, na.rm = TRUE),
    # LVAR
    SER_phi11_NFS = mean(NFS_phi11_se, na.rm = TRUE) / sd(NFS_phi11, na.rm = TRUE),
    SER_phi22_NFS = mean(NFS_phi22_se, na.rm = TRUE) / sd(NFS_phi22, na.rm = TRUE),
    SER_phi12_NFS = mean(NFS_phi12_se, na.rm = TRUE) / sd(NFS_phi12, na.rm = TRUE),
    SER_phi21_NFS = mean(NFS_phi21_se, na.rm = TRUE) / sd(NFS_phi21, na.rm = TRUE),
    SER_zeta1_NFS = mean(NFS_zeta1_se, na.rm = TRUE) / sd(NFS_zeta1, na.rm = TRUE),
    SER_zeta2_NFS = mean(NFS_zeta2_se, na.rm = TRUE) / sd(NFS_zeta2, na.rm = TRUE),
    SER_zeta12_NFS = mean(NFS_zeta12_se, na.rm = TRUE) / sd(NFS_zeta12, na.rm = TRUE),
    # NFS
    SER_phi11_SAM = mean(SAM_phi11_se, na.rm = TRUE) / sd(SAM_phi11, na.rm = TRUE),
    SER_phi22_SAM = mean(SAM_phi22_se, na.rm = TRUE) / sd(SAM_phi22, na.rm = TRUE),
    SER_phi12_SAM = mean(SAM_phi12_se, na.rm = TRUE) / sd(SAM_phi12, na.rm = TRUE),
    SER_phi21_SAM = mean(SAM_phi21_se, na.rm = TRUE) / sd(SAM_phi21, na.rm = TRUE),
    SER_zeta1_SAM = mean(SAM_zeta1_se, na.rm = TRUE) / sd(SAM_zeta1, na.rm = TRUE),
    SER_zeta2_SAM = mean(SAM_zeta2_se, na.rm = TRUE) / sd(SAM_zeta2, na.rm = TRUE),
    SER_zeta12_SAM = mean(SAM_zeta12_se, na.rm = TRUE) / sd(SAM_zeta12, na.rm = TRUE),
    # SAM
    SER_phi11_SEM = mean(SEM_phi11_se, na.rm = TRUE) / sd(SEM_phi11, na.rm = TRUE),
    SER_phi22_SEM = mean(SEM_phi22_se, na.rm = TRUE) / sd(SEM_phi22, na.rm = TRUE),
    SER_phi12_SEM = mean(SEM_phi12_se, na.rm = TRUE) / sd(SEM_phi12, na.rm = TRUE),
    SER_phi21_SEM = mean(SEM_phi21_se, na.rm = TRUE) / sd(SEM_phi21, na.rm = TRUE),
    SER_zeta1_SEM = mean(SEM_zeta1_se, na.rm = TRUE) / sd(SEM_zeta1, na.rm = TRUE),
    SER_zeta2_SEM = mean(SEM_zeta2_se, na.rm = TRUE) / sd(SEM_zeta2, na.rm = TRUE),
    SER_zeta12_SEM = mean(SEM_zeta12_se, na.rm = TRUE) / sd(SEM_zeta12, na.rm = TRUE),
    # SEM
    .groups = "drop")

#### performance overall ####
performance_overall <- performance %>% 
  summarise(across(AB_phi11_LVAR:SER_zeta12_SEM, ~ mean(.x, na.rm = TRUE)))

#### performance by effect size ####
performance_phi <- performance %>% 
  group_by(phi_size) %>% 
  summarise(across(AB_phi11_LVAR:SER_zeta12_SEM, ~ mean(.x, na.rm = TRUE)))

#### performance by n ####
performance_n <- performance %>% 
  group_by(n) %>% 
  summarise(across(AB_phi11_LVAR:SER_zeta12_SEM, ~ mean(.x, na.rm = TRUE)))

#### performance by obs ####
performance_obs <- performance %>% 
  group_by(obs) %>% 
  summarise(across(AB_phi11_LVAR:SER_zeta12_SEM, ~ mean(.x, na.rm = TRUE)))

#### performance by rho ####
performance_rho <- performance %>% 
  group_by(rho_gen) %>% 
  summarise(across(AB_phi11_LVAR:SER_zeta12_SEM, ~ mean(.x, na.rm = TRUE)))

