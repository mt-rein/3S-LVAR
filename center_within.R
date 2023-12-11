center_within <- function(data, vars, id){
  # data:
  #   a data frame
  # vars:
  #   a vector of variable names that will be within-person centered
  # id:
  #   a character that indicates the id variable (the variable that indicates
  #   which observations belong to which person)
  library(dplyr)
  data <- data %>% 
    group_by(!!rlang::sym(id)) %>% 
    mutate(across(all_of(vars), ~.x - mean(.x, na.rm = TRUE))) %>% 
    ungroup()
  
  return(data)
}