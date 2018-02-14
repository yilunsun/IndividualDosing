#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "Model.h"
#include "NodeTree.h"
#include "DeltaGrow.h"


DeltaGrow::DeltaGrow(const Node *oldLeaf,int newVar,double newTau) : Delta(DELTA_GROW)
{
  this->oldLeaf = oldLeaf;
  this->newVar = newVar;
  this->newTau = newTau;

  return;
}



DeltaGrow::~DeltaGrow(void)
{
  return;
};



NodeTree *DeltaGrow::PerformChange(NodeTree *tree) const
{
  Node *pp = (Node *) oldLeaf;

  Node *NewLeft = new Node;
  NewLeft->SetObservation(pp->GetObservation());
  NewLeft->SetSplitVariable(-1);
  NewLeft->SetMinimumLeafSize(pp->GetMinimumLeafSize());

  Node *NewRight = new Node;
  NewRight->SetObservation(pp->GetObservation());
  NewRight->SetSplitVariable(-1);
  NewRight->SetMinimumLeafSize(pp->GetMinimumLeafSize());

  pp->SetLeftNode(NewLeft);  pp->SetRightNode(NewRight);
  pp->SetSplitVariable(newVar);
  pp->SetSplitLevel(newTau);
  if (!pp->UpdateSubjectList())
    Rcout << "UpdateSubjectList failed in Delta Grow!"<<endl;    // this should never happen!

  return tree;
};
