#ifndef DELTABASIC_H
#define DELTABASIC_H

#ifndef _DIAGRAMS_DELTABASIC_H
#define _DIAGRAMS_DELTABASIC_H


#include "Delta.h"

class DeltaBasic : public Delta
{
 public:
  DeltaBasic(NodeTree *tree);
  ~DeltaBasic(void);

  NodeTree *PerformChange(NodeTree *tree) const;
  const NodeTree *GetNewTree(void) const {return newTree;};

 protected:
  NodeTree *newTree;
};


#endif
#endif