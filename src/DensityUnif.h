#ifndef DENSITYUNIF_H
#define DENSITYUNIF_H

#ifndef _DIAGRAMS_DENSITYUNIF_H
#define _DIAGRAMS_DENSITYUNIF_H

#include "Density.h"

class DensityUnif : public Density
{
public:
  DensityUnif(double Lower, double Upper);
  ~DensityUnif(void);

  double Sample(Random &ran) const;
  double Potential(double x,Random &ran) const;
  std::vector <double> PosteriorMean( NumericVector &data) const {std::vector <double> t; t.push_back(0); return t;};
  virtual std::vector <double> Predict(NumericVector &data) const {std::vector <double> t; t.push_back(0); return t;};

  Density *Copy(void) const;

private:
  double lower, upper;
};


#endif
#endif
