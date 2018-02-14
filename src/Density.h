#ifndef DENSITY_H
#define DENSITY_H

#ifndef _DIAGRAMS_DENSITY_H
#define _DIAGRAMS_DENSITY_H

#include <vector>
#include "Random.h"
#include <Rcpp.h>
using namespace Rcpp;

class Density
{
 public:
  Density(void);
  virtual ~Density(void);

  virtual double Sample(Random &ran) const = 0;
  virtual double Potential(double x, Random &ran) const = 0;

  virtual Density *Copy(void) const = 0;
  virtual std::vector <double> PosteriorMean(NumericVector &data) const=0;

  virtual std::vector <double> Predict(NumericVector &data) const=0;

 private:

};


inline Density::Density(void)
{
  return;
};


inline Density::~Density(void)
{
  return;
};


#endif
#endif
