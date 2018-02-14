#ifndef DELTACHANGE_H
#define DELTACHANGE_H

#ifndef _DIAGRAMS_DELTACHANGE_H
#define _DIAGRAMS_DELTACHANGE_H

#include "Delta.h"

class DeltaChange : public Delta
{
 public:
  DeltaChange(NodeTree *tree,const Node *nodeOld,const Node *nodeNew,int varOld,int varNew,double tauOld,double tauNew);
  ~DeltaChange(void);

  NodeTree *PerformChange(NodeTree *tree) const;
  const NodeTree *GetNewTree(void) const {return newTree;};

  const Node *GetNodeOld(void) const {return nodeOld;};
  const Node *GetNodeNew(void) const {return nodeNew;};

  int GetVarOld(void) const {return varOld;};
  int GetVarNew(void) const {return varNew;};
  double GetTauOld(void) const {return tauOld;};
  double GetTauNew(void) const {return tauNew;};
  
 protected:
  NodeTree *newTree;
  const Node *nodeOld;
  const Node *nodeNew;

  int varOld;
  int varNew;
  double tauOld;
  double tauNew;
};


#endif
#endif