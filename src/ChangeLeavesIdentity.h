#ifndef CHANGELEAVESIDENTITY_H
#define CHANGELEAVESIDENTITY_H

#ifndef _DIAGRAMS_CHANGELEAVESIDENTITY_H
#define _DIAGRAMS_CHANGELEAVESIDENTITY_H

#include "ChangeLeaves.h"

class ChangeLeavesIdentity : public ChangeLeaves
{
 public:
  ChangeLeavesIdentity(void);
  ~ChangeLeavesIdentity(void);

  double ProposeChange(std::vector<std::vector <int> > &leaves,const Observation *obs,
		       int MinimumLeafSize,Random &ran) {return 0.0;};
  double ProposeReverse(std::vector<std::vector <int> > &leaves,const Observation *obs,
			int MinimumLeafSize,Random &ran) {return 0.0;};
  ChangeLeaves *copy(void) const;

 protected:

};


inline ChangeLeavesIdentity::ChangeLeavesIdentity(void) : ChangeLeaves()
{
  return;
};


inline ChangeLeavesIdentity::~ChangeLeavesIdentity(void)
{
  return;
};



inline ChangeLeaves *ChangeLeavesIdentity::copy(void) const
{
  ChangeLeaves *r = new ChangeLeavesIdentity;

  return r;
};

#endif
#endif
