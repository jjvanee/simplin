#include <Rcpp.h>
using namespace Rcpp;

//' Do values in a numeric vector fall in specified range?
//'
//' This is a shortcut for `x >= left & x <= right`, implemented
//' efficiently in C++ for local values, and translated to the
//' appropriate SQL for remote tables.
//'
//' @param x A numeric vector of values
//' @param left,right Boundary values
//' @name simp_lin
//' @title Fast Simple Linear Regression`
//' @export
//' @examples
//' between(1:12, 7, 9)
//'
//' x <- rnorm(1e2)
//' x[between(x, -1, 1)]

// Create simp_lin function
// Anticipated output and inputs are List and Numeric Vectors respectively
// [[Rcpp::export]]
List simp_lin( NumericVector x, NumericVector y){
  // Extract size of dataset
  int n = x.size();
  // Allocate space for regression coefficients, vector length 2 with zeros 
  NumericVector coef (2);
  // Calculate intermediates 
  double var_x = var(x);
  double x_bar = mean(x);
  double cov_xy = 0;
  for(int i = 0; i < n; i++) {
    cov_xy += (x[i]-x_bar)*(y[i]-mean(y))/(n-1);
  }
  // Calculate Coefficients
  coef[1] = cov_xy/var_x;
  coef[0] = mean(y) - coef[1]*x_bar;
  // Checked code above against lm, it's good
  // Storage for Prediction Vector 
  NumericVector pred (n); 
  // Fill Prediction Vector
  for(int i = 0; i < n; i++) {
    pred[i] = coef[0]+coef[1]*x[i];
  }
  // Storage for residuals
  NumericVector resd (n);
  // Get residuals
  for(int i = 0; i < n; i++) {
    resd[i] = pred[i]-y[i];
  }
  // Derive MSE
  double MSE = 0;
  for(int i = 0; i < n; i++) {
    MSE += pow(resd[i], 2)/(n-2);
  }
  // Create Storage for Standard Errors
  NumericVector stde (2);
  stde[1] = sqrt(MSE/((n-1)*var_x));
  stde[0] = sqrt(MSE/n + pow(x_bar*stde[1], 2));
  // Quantiles of t-distribution 
  double quant = R::qt( 0.975, n-2, false, false);
  // Store CIs in Matrix
  NumericMatrix coef_CI(2,2);
  // Give row and column names to matrix
  //CharacterVector col = {"Intercept","Slope"};
  //CharacterVector row = {"0.975","0.025"};
  //colnames(coef_CI) = CharacterVector({"I","S"});
  //rownames(coef_CI) = CharacterVector({"U","L"});
  // Intercept
  coef_CI(0,0) = coef[0] + quant*stde[0];
  coef_CI(1,0) = coef[0] - quant*stde[0];
  // Slope Coefficient 
  coef_CI(0,1) = coef[1] + quant*stde[1];
  coef_CI(1,1) = coef[1] - quant*stde[1];
  // List with outputs
  List lin_out = List::create(coef, stde, coef_CI, pred, resd);
  return lin_out;
}

