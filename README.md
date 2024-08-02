This is a collection of functions that implement the method [Three-Step Latent Vector Autoregression (Rein, Vermunt, De Roover, Vogelsmeier, 2023)](https://osf.io/preprints/psyarxiv/a2muk) in R.

# Installation
The script `installation.R` provides some code to install or load the functions. See the file `Instructions.docx` for more information on how to load and use the functions. The related publication also provides an application example of the method using these functions.

# Required Packages
The functions are built on the `lavaan`, `dplyr`, `tidyr`, `RcppAlgos`, and `numDeriv` packages.

# Example Data
The `examples` folder contains some (simulated) data that can be used to test the functions.

# Main Functions
`step1(data, measurementmodel, id)`  
This function performs step 1 of 3S-LVAR.  
Input:
-	`data`: a data frame with an id variable and the indicator variables for the measurementmodel
-	`measurementmodel`: a string describing the measurement model using lavaan syntax
-	`id`: a character that indicates the id variable in data (the variable that indicates which observations belong to which person)

Output: a list with two elements
-	`fit`: the lavaan fit object
-	`data`: the input data frame (which is needed again for `step2()`)

`step2(step1output)`  
This function performs step 2 of 3S-LVAR.  
Input:
-	`step1output`: the object that was generated using the step1() function

Output: a list with four elements:
-	`data`: the input data frame to which the factor scores have been appended
-	`rho`: a vector containing the values of ρ (i.e., the model-based reliabilities)
-	`kappa`: a vector containing the values of κ
-	`fit_step1`: the lavaan fit object from step 1

`step3(step2output, structuralmodel = NULL)`  
This function performs step 3 of 3S-LVAR.  
Input:
-	`step2output`: the object that was generated using the step2() function
-	`structuralmodel`: optional, a string describing the structural model using lavaan syntax. Note: The lagged variables are created automatically by the function. They are named by appending “_lag” to the names of the latent constructs in step 1. For example, if the factors in step1() have been named “f1” and “f2”, then the lagged variables of the factor scores are automatically named “f1_lag” and “f2_lag”. These names must be adhered to when specifying the structuralmodel or the function will not work.

Output: a list with two elements:
-	`fit_step3`: the lavaan fit object
-	`data`: the data set that was used to estimate the model

`stepwiseSE(step2output, step3output)`  
This function adjusts the standard errors for the stepwise estimation as described in the article.  
Input:
-	`step2output`: the object that was generated using the step2() function
-	`step3output`: the object that was generated using the step3() function

Output: a list with three elements:
-	`SE`: the vector of adjusted standard errors
-	`z_values`: the vector of adjusted z-values
-	`p_values`: the vector of adjusted p-values

# Helper Functions
`center_within(data, vars, id)`  
This function centers the variables per person on their observed within-person mean.  
Input:
-	data: a data frame
-	vars: a vector of variable names that will be within-person centered
-	id: a character that indicates the id variable in data (the variable that indicates which observations belong to which person)

Output: the data frame in which the variables in `vars` have been within-person centered. Note: the variables have been overwritten (i.e., no new variables are created).

`correct_nickellsbias`  
This function adjusts the AR parameter(s) for Nickell's bias, if observed person-mean centered scores have been used, as described in the article.  
Input:
- `phi`: Either a vector with 1 element or a matrix with the AR parameters on the diagonal. The `phi` element of the output of `step3` can be used here.

Output:
- `phi`: The input after correcting for Nickell's bias.

