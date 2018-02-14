#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "NodeTree.h"
#include "Delta.h"
#include "DP.h"
#include "ModelSplitVariableIndependent.h"


ModelSplitVariableIndependent::ModelSplitVariableIndependent(const std::vector<double> &probability,
															 const std::vector<Density *> &density)
{
	prob = probability;

	int k;
	for (k = 0; k < density.size(); k++)
		dens.push_back(density[k]->Copy());

	return;
};



ModelSplitVariableIndependent::~ModelSplitVariableIndependent(void)
{
	int k;
	for (k = 0; k < dens.size(); k++)
		delete dens[k];

	return;
};



double ModelSplitVariableIndependent::PotentialDifference(const NodeTree &oldTree,
														  const Delta &newTree,
														  Random &ran) const
{
	double pot = DP::PotentialDifference(oldTree,newTree,*this,ran);

	return pot;
};





void ModelSplitVariableIndependent::Simulate(NodeTree &tree,Random &ran) const
{
  std::vector<Node *> node;

	node.push_back(&tree);
	int i;
	for (i = 0; i < node.size(); i++)
	{
		if (node[i]->GetRightNode() != NULL)
		{
			node.push_back(node[i]->GetRightNode());
			node.push_back(node[i]->GetLeftNode());

			int nr = ran.Discrete(prob);
			std::vector <int> SplittingVar = node[i]->GetSplittingVar();

			while (nr < node[i]->GetObservation()->GetNumCats()){//nr is dummy variable
        if (std::find(SplittingVar.begin(), SplittingVar.end(), nr) != SplittingVar.end()) {
          node[i]->GetRightNode()->ExcSplittingVar(nr);
          node[i]->GetLeftNode()->ExcSplittingVar(nr);
          break;
        } else {
          nr = ran.Discrete(prob);
        };
			};
			double tau = dens[nr]->Sample(ran);

			node[i]->SetSplitVariable(nr);
			node[i]->SetSplitLevel(tau);
		};
	};

	return;
};



Model *ModelSplitVariableIndependent::Copy(void) const
{
	Model *model = new ModelSplitVariableIndependent(prob,dens);

	return model;
};



double ModelSplitVariableIndependent::Potential(const NodeTree &tree, Random &ran) const
{
	double pot = 0.0;

  std::vector<const Node *> node;

	node.push_back(&tree);
	int i;
	for (i = 0; i < node.size(); i++)
	{
		if (node[i]->GetRightNode() != NULL)  // interior node
		{
			node.push_back(node[i]->GetRightNode());
			node.push_back(node[i]->GetLeftNode());

			int nr = node[i]->GetSplitVariable();
			double tau = node[i]->GetSplitLevel();
			pot += PotentialOneSplit(nr,tau,ran);
		}
	}

	return pot;
};



double ModelSplitVariableIndependent::PotentialOneSplit(int nr,double tau,Random &ran) const
{
	double pot = 0.0;

	pot += - log(prob[nr]);
	pot += dens[nr]->Potential(tau,ran);

	return pot;
};
