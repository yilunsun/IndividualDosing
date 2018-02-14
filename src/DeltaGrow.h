#ifndef DELTAGROW_H
#define DELTAGROW_H

#ifndef _DIAGRAMS_DELTAGROW_H
#define _DIAGRAMS_DELTAGROW_H

#include "Delta.h"

class DeltaGrow : public Delta
{
 public:
  DeltaGrow(const Node *oldLeaf,int newVar,double newTau);
  ~DeltaGrow(void);

  NodeTree *PerformChange(NodeTree *tree) const;
  const Node *GetOldLeaf(void) const {return oldLeaf;};
  int GetNewVar(void) const {return newVar;};
  double GetNewTau(void) const {return newTau;};

 protected:
  const Node *oldLeaf;
  int newVar;
  double newTau;
};


#endif
#endif