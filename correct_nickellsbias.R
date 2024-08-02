correct_nickellsbias <- function(phi, obs){
  if(is.vector(phi)){
    phi <- (phi*(obs-1) + 1)/(obs-2)
  }
  if(is.matrix(phi)){
    diag(phi) <- (diag(phi)*(obs-1) + 1)/(obs-2)
  }
  
  return(phi)
}