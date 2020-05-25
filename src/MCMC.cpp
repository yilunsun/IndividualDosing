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
						Random &ran, double T0, int Length, int i) const
{
  double Tt = pow(T0, 1-((double) i)/((double) Length));
	nAccept.resize(proposal.size());
	for (int j = 0; j < nAccept.size(); j++) nAccept[j] = 0;

	for (int i = 0; i < numberOfIteration; i++)
	{
	  //Rcout<<"In "<<i<<"-th iteration: \n";
		for (int j = 0; j < proposal.size(); j++)
		{
		  //Rcout<<"In "<<j<<"-th proposal:" <<tree->GetMiniNodeSize()<<"\n";
			double pot = 0.0;
			Delta *delta = proposal[j]->ProposeChange(*tree,pot,ran);
			//Rcout<<"point 1: "<<tree->GetMiniNodeSize()<<"\n";
			if (delta!=NULL)
			{
			// 	pot += priorStructure->PotentialDifference(*tree,*delta,ran);
			//   pot += priorSplit->PotentialDifference(*tree,*delta,ran);
			// 	pot += like->PotentialDifference(*tree,*delta,ran);
				//Rcout<<"point 2: "<<tree->GetMiniNodeSize()<<"\n";
        // if (pot < 0) pot = 0;
        //if (pot > 4) pot = pot/50; 
        //Rcout<<"aplha: "<<exp(- pot)<<"\t";
        NodeTree *tree_temp=tree->CopyTree();//Rcout<<"3. tree_temp: "<<tree_temp->GetMiniNodeSize()<<"\n";Rcout<<"3. tree: "<<tree->GetMiniNodeSize()<<"\n";
			  tree = delta->PerformChange(tree);
			  
			  double oldval = tree_temp->GetVal()["value"];
			  double newval = tree->GetVal()["value"];
			  
				if (ran.Unif01() <= exp((newval - oldval)/Tt)) // accept
				{
				  //Rcout<<"successful! \n";
					nAccept[j]++;//Rcout<<"success #: "<<nAccept[j]<<"\n";
					
					//keep a copy of temporary
					//Rcout<<"4. tree_temp: "<<tree_temp->GetMiniNodeSize()<<"\n";Rcout<<"4. tree: "<<tree->GetMiniNodeSize()<<"\n";
					
					// decide whether to revert the change
					if (tree->GetMiniNodeSize() >= tree->GetMinLeaf()) {
					  //tree = delta->PerformChange(tree);//
					  //Rcout<<"Success!"<<endl;
					  delete tree_temp;
					} else {
					  delete tree;
					  tree = tree_temp->CopyTree();
					};
					
					//delete tree_temp;
				};
				delete delta;
				delta=NULL;
			};
		};
	};

	return tree;
};
