#include <cmath>
#include "Random.h"
#include "DensityUnif.h"

DensityUnif::DensityUnif(double Lower, double Upper) : Density()
{
	lower=Lower;
	upper=Upper;
	return;
};

DensityUnif::~DensityUnif(void)
{
	return;
};

double DensityUnif::Sample(Random &ran) const
{
	double x = ran.Unif01();
	return x;
};

double DensityUnif::Potential(double x, Random &ran) const
{
	double pot = ran.PotentialUnif(x, lower, upper);

	return pot;
};

Density *DensityUnif::Copy(void) const
{
	Density *dens = new DensityUnif(lower, upper);

	return dens;
};
