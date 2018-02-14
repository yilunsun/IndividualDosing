#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "Model.h"
#include "NodeTree.h"
#include "DeltaChange.h"


DeltaChange::DeltaChange(NodeTree *tree,const Node *nodeOld,const Node *nodeNew,
			 int varOld,int varNew,double tauOld,double tauNew) : Delta(DELTA_CHANGE)
{
  newTree = tree;
  this->nodeOld = nodeOld;
  this->nodeNew = nodeNew;

  this->varOld = varOld;
  this->varNew = varNew;
  this->tauOld = tauOld;
  this->tauNew = tauNew;

  return;
}



DeltaChange::~DeltaChange(void)
{
  newTree->RemoveSubTree();

  delete newTree;
  
  return;
};



NodeTree *DeltaChange::PerformChange(NodeTree *tree) const
{
  tree->RemoveSubTree();
  delete tree;

  return newTree->CopyTree();
};
