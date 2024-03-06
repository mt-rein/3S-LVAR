# only for testing!
# data <- read.csv("examples/exampledata.csv")
# measurementmodel <- "
#   f1 =~ v1 + v2 + v3 + v4
#   f2 =~ v5 + v6 + v7 + v8
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
