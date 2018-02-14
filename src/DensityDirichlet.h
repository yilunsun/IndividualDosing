#ifndef DENSITYDIRICHLET_H
#define DENSITYDIRICHLET_H

#ifndef _DIAGRAMS_DENSITYDIRICHLET_H
#define _DIAGRAMS_DENSITYDIRICHLET_H

#include "Density.h"

class DensityDirichlet : public Density
{
public:
    DensityDirichlet(vector<double> alpha, int k);
    DensityDirichlet(void);
    ~DensityDirichlet(void);

    double Sample(Random &ran) const;
    double Potential(double x,Random &ran) const;
    Density *Copy(void) const;
    std::vector <double> PosteriorMean(NumericVector &data) const;
    std::vector <double> Predict(NumericVector &data) const;


private:
    vector<double> alpha;
    int k;

};


#endif
#endif
