#ifndef DELTAPRUNE_H
#define DELTAPRUNE_H

#ifndef _DIAGRAMS_DELTAPRUNE_H
#define _DIAGRAMS_DELTAPRUNE_H


#include "Delta.h"


class DeltaPrune : public Delta
{
 public:
  DeltaPrune(const Node *newLeaf);
  ~DeltaPrune(void);

  NodeTree *PerformChange(NodeTree *tree) const;
  const Node *GetNewLeaf(void) const {return newLeaf;};

 protected:
  const Node *newLeaf;
};


#endif
#endif