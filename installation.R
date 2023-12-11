# The following two lines install/update the packages required for the functions.
packages <- c("lavaan", "dplyr", "tidyr", "RcppAlgos", "numDeriv")
install.packages(packages)

# There are two ways to load the functions.
# 1) Download the function files from https://github.com/mt-rein/3S-LVAR and move 
# them into your working directory. Then, run the following lines:
source("step1.R")
source("step2.R")
source("step3.R")
source("stepwiseSE.R")
source("center_within.R")

# 2) Alternatively, you can load the functions directly from GitHub using the following
# lines of code:
install.packages("devtools") # required if you don't have the package installed yet
library(devtools)
#step1:
source_url("https://raw.githubusercontent.com/mt-rein/3S-LVAR/main/step1.R?token=GHSAT0AAAAAACLOXLBYE62NTPZRPNIYVQUKZLW76VA")
#step2:
source_url("https://raw.githubusercontent.com/mt-rein/3S-LVAR/main/step2.R?token=GHSAT0AAAAAACLOXLBY7I6E5HGULAY2SZV6ZLXAHCQ")
#step3:
source_url("https://raw.githubusercontent.com/mt-rein/3S-LVAR/main/step3.R?token=GHSAT0AAAAAACLOXLBZWRBNB7W553RSSQRYZLXAHPA")
#stepwiseSE:
source_url("https://raw.githubusercontent.com/mt-rein/3S-LVAR/main/stepwiseSE.R?token=GHSAT0AAAAAACLOXLBY5TJWIZYQKQRLC5OOZLXAINQ")
#center_within:
source_url("https://raw.githubusercontent.com/mt-rein/3S-LVAR/main/center_within.R?token=GHSAT0AAAAAACLOXLBZI7LBWGB6RFQASKKEZLXAJBQ")
