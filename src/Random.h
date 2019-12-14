#ifndef RANDOM_H
#define RANDOM_H

#ifndef _DIAGRAMS_RANDOM_H
#define _DIAGRAMS_RANDOM_H

#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;


#define MULTIPLIER 69069
#define SHIFT          1
#define MODULUS    256*256*256*128
#define INVMOD     ( (double) 1 / ((double) MODULUS)) / ((double) 2)


class Random
{
 public:
  Random(unsigned int seed);
  Random();
  ~Random(void);

  double Unif(void);
  double Unif01(void);
  double Norm01(void);
  double Exponential(double lambda);
  int Poisson(double lambda);
  int Binomial(int n,double p);
  int Discrete(const std::vector<double> &prob);
  double Gamma(int alpha, double beta);

  double PotentialGaussian(double variance,double mean,double x);
  double PotentialExp(double lambda,double x);
  double PotentialWeibull(double alpha,double beta, double x);
  double PotentialPoisson(double lambda,int x);
  double PotentialUnif(double x, double lower, double upper);
  double PotentialUnifDiscrete(double x, int lower, int upper);
  double PotentialTruncatedPoisson(double lambda,int x, int length);
  double PotentialBinomial(int n,double p,int x);
  double PotentialTruncatedGaussian(double variance,double mean,double x, double lower, double upper);
  
  double lnGamma(double x);
    unsigned int getseed(){return seedValue;};

 private:
  unsigned int seedValue;
  void setseed();
  int haveNorm01;
  double norm;


  double erfinv(double y);
  double erf(double x);
  double normcdf(double x, double mu, double sd);
  double fix(double x);
  double erfcore(double x);
  double abswu(double x);
//  vector <double> Random::normrand(int Size, double mu, double sd);
//  double Random::truncnorm(double mu, double sigma, double lower, double upper);
  double norminv(double p, double mu, double sd);

};

#endif
#endif
