# only for testing!
# data <- read.csv("exampledata.csv")
# measurementmodel <- "
#   f1 =~ V1 + V2 + V3 + V4 + V5
#   f2 =~ V6 + V7 + V8 + V9 + V10
# 
#   V1 ~ 0*1
#   V2 ~ 0*1
#   V3 ~ 0*1
#   V4 ~ 0*1
#   V5 ~ 0*1
#   V6 ~ 0*1
#   V7 ~ 0*1
#   V8 ~ 0*1
#   V9 ~ 0*1
#   V10 ~ 0*1
#   "
# id = "id"

step1 <- function(data, measurementmodel, id){
  # data:
  #   a data frame with the indicator and ID variables
  # measurementmodel:
  #   a string describing the measurement model using the lavaan syntax
  # id:
  #   a character that indicates the id variable (the variable that indicates
  #   which observations belong to which person)
  
  library(lavaan)
  
  # estimate the measurement model:
  fit <- cfa(measurementmodel,
             data = data,
             orthogonal = TRUE,
             missing = "ML",
             cluster = id)

  # assemble output
  output <- list("fit_step1" = fit,
                 "data" = data)
  return(output)
}
