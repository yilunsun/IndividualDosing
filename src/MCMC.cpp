#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "Model.h"
#include "Delta.h"
#include "NodeTree.h"
#include "Proposal.h"
#include "ProposalBasicRadical.h"
#include "MCMC.h"


MCMC::MCMC(Model &priorTreeStructure,Model &priorSplitVariable,Model &likelihood)
{
	priorStructure = priorTreeStructure.Copy();
	priorSplit = priorSplitVariable.Copy();
	like = likelihood.Copy();

	return;
};



MCMC::~MCMC(void)
{
	delete priorStructure;
	delete priorSplit;
	delete like;

	return;
};



NodeTree *MCMC::Iterate(NodeTree *tree,std::vector<Proposal *> proposal,
                        int numberOfIteration,std::vector<int> &nAccept,
                        Random &ran, double T0, const std::vector<double> &proprob, double complexity, bool verbose, bool np) const
{
  nAccept.resize(proposal.size());
  for (int j = 0; j < nAccept.size(); j++) nAccept[j] = 0;
  double value_best = -9999;
  NodeTree *tree_best;
  for (int i = 0; i < numberOfIteration; i++)
  {
    if (verbose)
      Rcout<<"Inner Loop: "<<i<<endl;
    
    double Tt = pow(T0, 1-((double) i)/((double) numberOfIteration));
    // Rcout<<"Tt: "<<Tt<<endl;
    
    int j = ran.Discrete(proprob);
    
    double pot = 0.0;
    Delta *delta = proposal[j]->ProposeChange(*tree,pot,ran);
    //Rcout<<"point 1: "<<tree->GetMiniNodeSize()<<"\n";
    if (delta!=NULL)
    {
      NodeTree *tree_temp=tree->CopyTree();//Rcout<<"3. tree_temp: "<<tree_temp->GetMiniNodeSize()<<"\n";Rcout<<"3. tree: "<<tree->GetMiniNodeSize()<<"\n";
      tree = delta->PerformChange(tree);
      
      if (tree->GetMiniNodeSize() >= tree->GetMinLeaf()) // size is ok
      {
        List oldlist = tree_temp->GetVal(np);
        List newlist = tree->GetVal(np);
        
        double oldval = as<NumericVector>(oldlist["value"])[0] - tree_temp->GetSize() * complexity;
        double newval = as<NumericVector>(newlist["value"])[0] - tree->GetSize() * complexity;
        
        if (ran.Unif01() <= exp((newval - oldval)/Tt)) // accept
        {
          if (newval > value_best)
          {
            value_best = newval;
            tree_best = tree->CopyTree();
          }
          nAccept[j]++;
          delete tree_temp;
        } else {
          delete tree;
          tree = tree_temp->CopyTree(); //revert move if accept is not OK
        }
      } else {
        delete tree;
        tree = tree_temp->CopyTree(); //revert move if size is not OK
      }
      
      delete delta;
      delta=NULL;
    };
  };
  
  if (value_best > tree->GetValue()) {
    delete tree;
    tree = tree_best->CopyTree();
    delete tree_best;
    tree_best = NULL;
  } else {
    delete tree_best;
    tree_best = NULL;
  }
  
  return tree;
};
