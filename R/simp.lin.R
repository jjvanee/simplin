library(Rcpp)
simp_lin_R <- function(x, y){
  ### Check if vectors are if same length
  if(length(x) != length(y)){stop( "x and y are of differing lengths!" )}
  ### Check if vectors are numeric
  if(!is.numeric(x) | !is.numeric(y)){stop( "x or y are not numeric!" )}
  simp_lin(x,y)
}
