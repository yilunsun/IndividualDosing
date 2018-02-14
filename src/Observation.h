#ifndef OBSERVATION_H
#define OBSERVATION_H

#ifndef _DIAGRAMS_OBSERVATION_H
#define _DIAGRAMS_OBSERVATION_H

#include <vector>
#include <string>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

class Observation
{
 public:
  Observation(const NumericMatrix &x, const IntegerVector &y, const int &cat_num, const NumericMatrix &V);
  Observation();
  ~Observation(void);

  int GetP(void) const {return p;};
  int GetN(void) const {return n;};
  int GetK(void) const {return k;};
  int GetY(int i) const {return y[i];};
  double GetX(int i,int j) const {return x(i,j);};
  double GetV(int i,int j) const {return V(i,j);};
  int GetNumCats(void) const {return cat_num;}; //number of categorical variables
  void Show();
  void SummaryStat();

 private:
  int p; // number of predictors
  int n; // number of observations
  int k; // number of categories in outcome

  IntegerVector y;
  NumericMatrix x;
  int cat_num;
  NumericMatrix V;// value matrix
};

#endif
#endif
