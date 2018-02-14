#include "Random.h"
#include "Observation.h"
#include "Model.h"
#include "Density.h"
#include "NodeTree.h"
#include "Delta.h"
#include "DeltaChange.h"
#include "ProposalChange.h"


ProposalChange::ProposalChange(const std::vector<double> &probability,
			       const std::vector<Density *> &density)
{
  prob = probability;

  int k;
  for (k = 0; k < density.size(); k++)
    dens.push_back(density[k]->Copy());

  return;
};


ProposalChange::~ProposalChange(void)
{
  for (int k = 0; k < dens.size(); k++)
    delete dens[k];

  return;
};


Delta *ProposalChange::ProposeChange(const NodeTree &tree,double &potential,Random &ran) const
{
  NodeTree *newTree = tree.CopyTree();

  if (tree.GetRightNode() == NULL)    // no interior nodes exist!
    {                                 // propose an unchanged tree.
      Delta *delta = NULL;
      potential = 0.0;                // no interior node also for reverse move.
      newTree->RemoveSubTree();
	    delete newTree;
	    newTree=NULL;
      return delta;
    };


  //
  // Now we know the root is interior. Find all interior nodes
  //

  std::vector<Node *> node;
  std::vector<const Node *> nodeOld;
  node.push_back(newTree);
  nodeOld.push_back(&tree);
  int i;
  for (i = 0; i < node.size(); i++)
  {
    if (node[i]->GetRightNode()->GetLeftNode() != NULL)
      {
        node.push_back(node[i]->GetRightNode());
        nodeOld.push_back(nodeOld[i]->GetRightNode());
      }

    if (node[i]->GetLeftNode()->GetRightNode() != NULL)
      {
        node.push_back(node[i]->GetLeftNode());
        nodeOld.push_back(nodeOld[i]->GetLeftNode());
      }
  };

  //
  // draw one of the interior nodes uniformly
  //
  std::vector<double> pr = std::vector<double>(node.size(),1.0/node.size());
  int nr = ran.Discrete(pr);  // i.e. propose new rule for this node

  //
  // draw new variable number and tau
  //

  int newVar = ran.Discrete(prob);
  double newTau = dens[newVar]->Sample(ran);
  int oldVar = node[nr]->GetSplitVariable();
  double oldTau = node[nr]->GetSplitLevel();

  node[nr]->SetSplitVariable(newVar);
  node[nr]->SetSplitLevel(newTau);

  potential = 0.0;
  potential -= -log((double) prob[newVar]);
  potential -= dens[newVar]->Potential(newTau,ran);

  potential += -log((double) prob[oldVar]);
  potential += dens[oldVar]->Potential(oldTau,ran);

  //
  // Then we have to update the subject list for all
  // descendent of the chosen node. Let a function
  // in the Node do this.
  //

  if (!node[nr]->UpdateSubjectList())
  {
    Delta *delta=NULL;
    newTree->RemoveSubTree();
    delete newTree;
    newTree=NULL;
    return delta;
  }
  else
  {
    Delta *delta = new DeltaChange(newTree,nodeOld[nr],node[nr],oldVar,newVar,oldTau,newTau);
    return delta;
  };
};
