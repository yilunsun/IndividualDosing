#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]


arma::colvec repc(const NumericVector& input, const int& ct){
  
  arma::colvec input_a = as<arma::colvec>(input);
  
  arma::colvec result = input_a;
  
  for (int i=0;i<ct-1;i++){
    result = arma::join_cols(result, input_a);
  }
  
  return result;
}

NumericVector seq_length(const double& start, const double& end, const int& length_out){
  NumericVector result;
  double increment = (end-start)/(length_out - 1);
  for (int i = 0; i < length_out; i++){
    result.push_back(start+i*increment);
  }
  return result;
}

NumericVector dirtyGBM(const NumericMatrix& x, const NumericVector& y, const NumericMatrix& newdf){ // this implementation requires expose predict in xgboost package as predict2
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:xgboost"); 
  
  // Make function callable from C++
  Rcpp::Function xgbm = package_env["xgboost"];  
  Rcpp::Function predict2 = package_env["predict2"]; 
  
  // Call the function and receive output (might not be list)
  //Rcout << "here";
  SEXP rf_obj = xgbm(Named("data")=x, _["label"]=y, _["max.depth"] =3 , _["eta"] = 1, _["nthread"] = 2, _["nrounds"] = 2, _["verbose"] = 0);
  //Rcout << "2";
  NumericVector yhat = predict2(rf_obj, newdf);
  
  return yhat;
}

List dirtyDensity(const arma::colvec& x){ // this implementation requires expose predict in xgboost package as predict2
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:stats"); 
  
  // Make function callable from C++
  Rcpp::Function densi = package_env["density"];
  
  //Rcout << "here";
  List dens = densi(x);
  
  return dens;
}

arma::colvec dirtySpline(const arma::colvec& x, const arma::colvec& y, const arma::colvec& newdf){ // this implementation requires expose predict in xgboost package as predict2
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:stats"); 
  
  // Make function callable from C++
  Rcpp::Function smoothspline = package_env["smooth.spline"];  
  Rcpp::Function predict2 = package_env["predict"]; 
  
  // Call the function and receive output (might not be list)
  //Rcout << "here";
  SEXP rf_obj = smoothspline(Named("x")=x, _["y"]=y);
  //Rcout << "2";
  List yhat = predict2(rf_obj, newdf);
  
  return yhat["y"];
}

arma::colvec w_fn(const double& bw, const arma::colvec a_vals, const arma::colvec& a, const int& n) {
  arma::colvec w_avals(a_vals.n_elem);
  for (int i=0; i<a_vals.n_elem; i++) {
    arma::colvec a_std = (a - a_vals(i))/bw;
    arma::colvec kern_std = arma::normpdf(a_std)/bw;
    w_avals(i) = arma::mean(pow(a_std,2.0) % kern_std) * (arma::normpdf(0.0)/bw)/(arma::mean(kern_std) * arma::mean(pow(a_std,2.0) % kern_std) - pow(arma::mean(a_std % kern_std),2.0));
  }
  return w_avals/n;
}

arma::colvec dirtyLocpoly(const double& h, const arma::colvec& out, const arma::colvec& a, const arma::colvec& a_out){
  
  Rcpp::Environment package_env("package:KernSmooth"); 
  
  Rcpp::Function locpoly = package_env["locpoly"];
  // 
  List polybw = locpoly(Named("x")=a, _["y"]=out, _["bandwidth"] = h);
  
  arma::colvec result;
  
  arma::interp1(as<arma::colvec>(polybw["x"]), as<arma::colvec>(polybw["y"]), a_out, result, "linear");
  
  return result;
}

double risk_fn(const double& h, const arma::colvec& a_vals, const arma::colvec& out, const arma::colvec& a, const int& n) {
  
  arma::colvec hats;
  
  arma::colvec a_interp = a;
  a_interp.elem( find(a_interp > max(a_vals)) ).fill(max(a_vals));
  a_interp.elem( find(a_interp < min(a_vals)) ).fill(min(a_vals));
  
  arma::interp1(a_vals, w_fn(h, a_vals, a, n), a_interp, hats, "linear");
  
  arma::colvec cts_eff_fn = dirtyLocpoly(h, out, a, a);
  // 
  arma::colvec one1(hats.n_elem);
  one1.ones();
  
  return arma::mean(pow((out - cts_eff_fn)/(one1 - hats),2.0));
}


// [[Rcpp::export]]
void vector_max(std::vector<double> v, double &max, int &imax){
  std::vector<double>::size_type p=0;
  imax = -1;
  max = std::numeric_limits<double>::lowest();
  
  for (auto &val : v)
  {
    if (!std::isnan(val) && val>max)
    {
      imax = p;
      max = val;
    }
    p++;
  }
}


// [[Rcpp::export]]
std::vector<double> dose_kernel_rf_cpp(const NumericVector& y, const NumericVector& a, const NumericMatrix& x, const NumericVector& a_out,  
                                                 double a_min = -999, double a_max = -999, const int n_pts = 100) 
{
  
  //require(KernSmooth)
  
  //kern <- function(t) {
  //  dnorm(t)
  //}
  int n = x.nrow();
  
  NumericVector bw_seq = seq_length(0.2, 10, 100);
  
  if (a_min == -999) 
  {
    a_min = min(a);
  }
  if (a_max == -999) 
  {
    a_max = max(a);
  }
  
  NumericVector a_vals = seq_length(a_min, a_max, n_pts);
  
  NumericVector pimod_vals = dirtyGBM(x, a, x);
  
  NumericVector pi2mod_vals = dirtyGBM(x, Rcpp::pow(a - pimod_vals, 2.0), x);
  
  arma::mat x_p = as<arma::mat>(x);
  
  arma::mat x_new = x_p;
  
  for (int i=0; i<n_pts; i++){
    x_new = arma::join_cols(x_new, x_p);
  }
  
  arma::colvec a_new(a_vals.size()*n);
  for (int i=0;i<a_vals.size();i++){
    for (int j=0;j<n;j++){
      a_new(i*n+j) = a_vals(i);
    }
  }
  
  arma::mat xa_new = arma::join_rows(x_new, arma::join_cols(as<arma::colvec>(a), a_new));
  
  arma::colvec muhat_vals = as<arma::colvec>(dirtyGBM(Rcpp::cbind(x,a), y, wrap(xa_new)));
  //Rcout<<"here aaa1 \n";
  arma::colvec a_std = (xa_new.col(xa_new.n_cols-1) - repc(pimod_vals,n_pts+1))/sqrt(abs(repc(pi2mod_vals,n_pts+1)));
  //Rcout<<"here aaa2 \n";
  arma::colvec pihat_vals;
  
  List density = dirtyDensity(a_std);
  
  arma::interp1(as<arma::colvec>(density["x"]), as<arma::colvec>(density["y"]), a_std, pihat_vals, "linear");
  //Rcout<<"here aaa3 \n";
  arma::colvec pihat = pihat_vals.subvec(0,n-1);
  
  arma::mat pihat_mat = arma::reshape(pihat_vals.subvec(n, pihat_vals.n_elem-1), n, n_pts);
  
  arma::colvec varpihat = dirtySpline(as<arma::colvec>(a_vals), arma::conv_to<arma::colvec>::from(arma::mean(pihat_mat, 0)), as<arma::colvec>(a));
  //Rcout<<"here aaa4 \n";
  arma::mat varpihat_mat;
  
  for (int i=0;i<n;i++){
    varpihat_mat.insert_rows(i,arma::conv_to<arma::colvec>::from(arma::mean(pihat_mat, 0)));
  }
  
  arma::colvec muhat = muhat_vals.subvec(0,n-1);
  
  arma::mat muhat_mat = arma::reshape(muhat_vals.subvec(n, muhat_vals.n_elem-1), n, n_pts);
  
  arma::colvec mhat = dirtySpline(as<arma::colvec>(a_vals), arma::conv_to<arma::colvec>::from(arma::mean(muhat_mat, 0)), as<arma::colvec>(a));
  //Rcout<<"here aaa5 \n";
  arma::mat mhat_mat;
  
  for (int i=0;i<n;i++){
    mhat_mat.insert_rows(i, arma::conv_to<arma::colvec>::from(arma::mean(muhat_mat, 0)));
  }
  
  arma::colvec pseudo_out = (as<arma::colvec>(y) - muhat)/(pihat/varpihat) + mhat;
  //Rcout<<"here aaa6 \n";
  arma::colvec risk_est;
  risk_est.ones(bw_seq.size());
  //Rcout<<"size of bw_seq: "<<bw_seq.size()<<"\n";
  //Rcout<<"size of risk_est: "<<risk_est.size()<<"\n";
  for (int i=0;i<risk_est.size();i++){
    risk_est(i) = risk_fn(bw_seq(i), a_vals, pseudo_out, as<arma::colvec>(a), n);
  }  
  //Rcout<<"here aaa7 \n";
  double h_opt = bw_seq(index_min(risk_est));
  
  arma::colvec est = dirtyLocpoly(h_opt, pseudo_out, a, a_out);
  
  return arma::conv_to<std::vector<double> >::from(est);
}