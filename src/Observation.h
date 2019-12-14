#ifndef OBSERVATION_H
#define OBSERVATION_H

#ifndef _DIAGRAMS_OBSERVATION_H
#define _DIAGRAMS_OBSERVATION_H

#include <vector>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;

class Observation
{
 public:
  Observation(const NumericMatrix& xo, const NumericMatrix& x, const IntegerVector& y, const int& cat_num, const NumericVector& V, const NumericVector& a, const NumericVector& candidate_dose);
  Observation();
  ~Observation(void);

  int GetP(void) const {return p;};
  int GetN(void) const {return n;};
  int GetK(void) const {return k;};
  int GetY(int i) const {return y[i];};
  double GetX(int i,int j) const {return x(i,j);};
  NumericVector GetXi(int i) const {return xo(i,_);};
  double GetV(int i) const {return V[i];};
  double GetA(int i) const {return a[i];};
  NumericVector GetCanDose(void) const {return candidate_dose;};
  int GetNumCats(void) const {return cat_num;}; //number of categorical variables
  void Show();
  void SummaryStat();

 private:
  int p; // number of predictors
  int n; // number of observations
  int k; // number of categories in outcome

  IntegerVector y;
  NumericMatrix x;
  NumericMatrix xo;
  int cat_num;
  NumericVector V;// value vector
  NumericVector a;//observed dose
  NumericVector candidate_dose;

};

#endif
#endif
