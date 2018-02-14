#ifndef DENSITYDUMMY_H_
#define DENSITYDUMMY_H_

#ifndef _DIAGRAMS_DENSITYDUMMY_H_
#define _DIAGRAMS_DENSITYDUMMY_H_

#include "Density.h"

class DensityDummy : public Density
{
public:
    DensityDummy(void);
    ~DensityDummy(void);

    double Sample(Random &ran) const {return 0.5;};
    double Potential(double x,Random &ran) const {return 0;};
    std::vector <double> PosteriorMean( NumericVector &data) const {std::vector <double> t; t.push_back(0); return t;};
    virtual std::vector <double> Predict(NumericVector &data) const {std::vector <double> t; t.push_back(0); return t;};

    Density *Copy(void) const;

private:
};


#endif
#endif
