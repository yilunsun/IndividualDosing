#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "Model.h"
#include "NodeTree.h"
#include "DeltaPrune.h"


DeltaPrune::DeltaPrune(const Node *newLeaf) : Delta(DELTA_PRUNE)
{
  this->newLeaf = newLeaf;

  return;
}



DeltaPrune::~DeltaPrune(void)
{
  return;
};



NodeTree *DeltaPrune::PerformChange(NodeTree *tree) const
{
  Node *pp = (Node *) newLeaf;

  delete pp->GetRightNode();
  delete pp->GetLeftNode();
  pp->SetLeftNode(NULL);
  pp->SetRightNode(NULL);
  pp->SetSplitVariable(-1);

  return tree;
};
