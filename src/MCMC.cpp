#include <cmath>
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
						Random &ran) const
{

	int j;
	nAccept.resize(proposal.size());
	for (j = 0; j < nAccept.size(); j++) nAccept[j] = 0;

	int i;
	for (i = 0; i < numberOfIteration; i++)
	{
		int j;
		for (j = 0; j < proposal.size(); j++)
		{
			double pot = 0.0;
			Delta *delta = proposal[j]->ProposeChange(*tree,pot,ran);
			if (delta!=NULL)
			{
				pot += priorStructure->PotentialDifference(*tree,*delta,ran);

				pot += priorSplit->PotentialDifference(*tree,*delta,ran);

				pot += like->PotentialDifference(*tree,*delta,ran);

				if (ran.Unif01() <= exp(- pot))
				{
					nAccept[j]++;
					tree = delta->PerformChange(tree);
					// For debugging only
					// Dump trees
					//string tempFilename("../temptree.txt");
					//tree->DumpDotFile(tempFilename, 0);
				}
				delete delta;
				delta=NULL;
			};
		}
	}

	return tree;
};
