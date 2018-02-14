#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "Model.h"
#include "NodeTree.h"
#include "Delta.h"
#include "DeltaBasic.h"
#include "ProposalSwap.h"


ProposalSwap::ProposalSwap()
{
	return;
};

ProposalSwap::~ProposalSwap(void)
{
	return;
};


Delta *ProposalSwap::ProposeChange(const NodeTree &tree,double &potential,
										Random &ran) const
{

	// Single root tree, no possible swap can be made
	// Return.
	if (tree.GetRightNode() == NULL)
	{
		Delta *delta=NULL;
		return delta;
	};

	// If the tree is of size two, no possible swap move can be made.
	// Return.

	if (tree.GetRightNode()->GetRightNode() == NULL && tree.GetLeftNode()->GetRightNode() == NULL)
	{
		Delta *delta=NULL;
		return delta;
	};


	//
	// Find all possible swappable pair.
	//
	NodeTree* NewTree=tree.CopyTree();

	// To grow the tree
	//
	// draw one of the leaf nodes uniformly
	//
	std::vector<Node *> node;
	std::vector<Node *> leaves;

	// PairA store one part of the possible pair
	// PairB store the corresponding the other part.
	std::vector<Node *> PairA;
	std::vector<Node *> PairB;
	//double potential = 0;
	node.push_back(NewTree);
	int i;
	for (i = 0; i < node.size(); i++)
	{
		if (node[i]->GetRightNode()==NULL)
			leaves.push_back(node[i]);
		else
		{
			if (node[i]->GetRightNode()->GetRightNode() != NULL)
			{
				PairA.push_back(node[i]);
				PairB.push_back(node[i]->GetRightNode());
			};
			if (node[i]->GetLeftNode()->GetRightNode() != NULL)
			{
				PairA.push_back(node[i]);
				PairB.push_back(node[i]->GetLeftNode());
			};
			node.push_back(node[i]->GetRightNode());
			node.push_back(node[i]->GetLeftNode());
		};


	}



	std::vector<double> pr = std::vector<double>(PairA.size(),1.0/PairA.size());
	int nr = ran.Discrete(pr);  // Among all the possible pairs, uniformly pick one.


	//Swap the splitting rules.

	int TempVar = PairA[nr]->GetSplitVariable();
	double TempLevel = PairA[nr]->GetSplitLevel();

	PairA[nr]->SetSplitVariable(PairB[nr]->GetSplitVariable());
	PairA[nr]->SetSplitLevel(PairB[nr]->GetSplitLevel());

	PairB[nr]->SetSplitVariable(TempVar);
	PairB[nr]->SetSplitLevel(TempLevel);


	potential = 0.0;

	//
	// Then we have to update the subject list for all
	// descendent of the chosen node. Let a function
	// in the Node do this.
	//

	if (!PairA[nr]->UpdateSubjectList())
	{
		Delta *delta=NULL;
		NewTree->RemoveSubTree();
		delete NewTree;
		NewTree=NULL;

		return delta;
	}
	else
	{
		Delta *delta = new DeltaBasic(NewTree);
		//NewTree->RemoveSubTree();
		//delete NewTree;
		//NewTree = NULL;

		return delta;
	};


	Delta *delta = NULL;
	potential = 0.0;

	return delta;
};
