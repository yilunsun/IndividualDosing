#include <cmath>
#include <numeric>
#include "Random.h"
#include "DensityDirichlet.h"

DensityDirichlet::DensityDirichlet(std::vector<double> alpha, int k) : Density()
{
    this->alpha=alpha;
    this->k=k;
    return;
};

DensityDirichlet::DensityDirichlet() : Density()
{
    return;
};

DensityDirichlet::~DensityDirichlet(void)
{
    return;
};

double DensityDirichlet::Sample(Random &ran) const
{
    double x = 0.0; // To be coded;
    return x;
};



double DensityDirichlet::Potential(double x,Random &ran) const
{
    double pot = 0.0; // To be coded;

    return pot;
};

Density *DensityDirichlet::Copy(void) const
{
    Density *dens = new DensityDirichlet(alpha,k);

    return dens;
};

std::vector <double> DensityDirichlet::PosteriorMean(NumericVector &data) const
{
  int n = data.size();
  int k = alpha.size();
  std::vector <double> Result;
  IntegerVector data1 = as<IntegerVector>(data);
  double alpha_sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);

  for(int i=1;i<k;i++)
  {
    int count = 0;

    double s1 = 0.0;

    for (int s = 0; s < n; s++)
    {
      if (data1[s]==i) count++;
    }

    s1 = (double(count)+alpha[i])/(alpha_sum+double(n));

    Result.push_back(s1);
  }

  return Result;
};


std::vector <double> DensityDirichlet::Predict(NumericVector &data) const
{
  std::vector <double> Result=PosteriorMean(data);

    return Result;
};
