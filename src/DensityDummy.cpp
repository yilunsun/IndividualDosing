#include <cmath>
#include "Random.h"
#include "DensityDummy.h"

DensityDummy::DensityDummy(void) : Density()
{
  return;
};


DensityDummy::~DensityDummy(void)
{
  return;
};


Density *DensityDummy::Copy(void) const
{
  Density *dens = new DensityDummy();

  return dens;
};
