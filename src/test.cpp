#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

/*** R
npCBPS_neo <- function(y,x) {
  library(CBPS)
  #print("X:", dim(X))
  return(1/(length(y)*unlist(npCBPS(y~x)['weights'])))
}

gam_neo <- function(a,ps,y,a_out) {
  library(mgcv)
  fit <- gam(y~s(a)+s(ps,by=a))
  pred <- predict(fit, data.frame(ps=ps,a=rep(a_out,each=length(ps))))
  optvalue <- apply(matrix(pred, nrow = length(ps)), 2, mean)
  return(optvalue)
}
*/
// [[Rcpp::export]]
NumericVector GCBPS(const NumericMatrix& x, const NumericVector& y){ // this implementation requires expose predict in xgboost package as predict2
  
  // Obtain environment containing function
  // Rcpp::Environment package_env("package:bypassFormula"); 
  // Rcpp::Function npCBPS_neo = package_env["npCBPS_neo"];  
  // 
  Rcpp::Environment G = Rcpp::Environment::global_env();
  Rcpp::Function npCBPS_neo = G["npCBPS_neo"];
  // Call the function and receive output (might not be list)
  //Rcout << "here";
  NumericVector rf_obj = npCBPS_neo(Named("y")=y, _["x"]=x);
  //Rcout << "2";
  
  return rf_obj;
}

// [[Rcpp::export]]
NumericVector GAM(const NumericVector& a, const NumericVector& ps, const NumericVector& y, const NumericVector& a_out){ 
  
  // Obtain environment containing function
  // Rcpp::Environment package_env("package:bypassFormula"); 
  // Rcpp::Function gam_neo = package_env["gam_neo"];  
  // 
  Rcpp::Environment G = Rcpp::Environment::global_env();
  Rcpp::Function gam_neo = G["gam_neo"];
  
  // Call the function and receive output (might not be list)
  //Rcout << "here";
  NumericVector a_pred = gam_neo(Named("a")=a, _["ps"]=ps, _["y"]=y, _["a_out"]=a_out);
  //Rcout << "2";
  
  return a_pred;
}
