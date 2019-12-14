#ifndef CHANGELEAVES_H
#define CHANGELEAVES_H

#ifndef _DIAGRAMS_CHANGELEAVES_H
#define _DIAGRAMS_CHANGELEAVES_H

#include "Observation.h"
#include "Random.h"

class ChangeLeaves
{
 public:
  ChangeLeaves(void);
  virtual ~ChangeLeaves(void);

  virtual double ProposeChange(std::vector<std::vector <int> > &leaves,const Observation *obs,
			       int MinimumLeafSize,Random &ran) = 0;
  virtual double ProposeReverse(std::vector<std::vector <int> > &leaves,const Observation *obs,
				int MinimumLeafSize,Random &ran) = 0;
  virtual ChangeLeaves *copy(void) const = 0;

 private:
};


inline ChangeLeaves::ChangeLeaves(void)
{
  return;
};


inline ChangeLeaves::~ChangeLeaves(void)
{
  return;
};


#endif
#endif
