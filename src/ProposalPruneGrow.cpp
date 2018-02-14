#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "Model.h"
#include "NodeTree.h"
#include "Delta.h"
#include "DeltaPrune.h"
#include "DeltaGrow.h"

#include "ProposalPruneGrow.h"



ProposalPruneGrow::ProposalPruneGrow(const std::vector<double> &probability,
									 const std::vector<Density *> &density)
{
	prob = probability;

	int k;
	for (k = 0; k < density.size(); k++)
		dens.push_back(density[k]->Copy());

	return;
};

ProposalPruneGrow::~ProposalPruneGrow(void)
{
	return;
};


Delta *ProposalPruneGrow::ProposeChange(const NodeTree &tree,double &potential,
										Random &ran) const
{
	if (tree.GetRightNode() == NULL)    // no interior nodes exist!
	{                                 // Grow the tree
		NodeTree *newTree = tree.CopyTree();
		Node *NewLeft=new Node;
		NewLeft->SetObservation(newTree->GetObservation());
		Node *NewRight=new Node;
		NewRight->SetObservation(newTree->GetObservation());
		NewLeft->SetSplitVariable(-1);
		NewRight->SetSplitVariable(-1);
		NewLeft->SetMinimumLeafSize(newTree->GetMinimumLeafSize());
		NewRight->SetMinimumLeafSize(newTree->GetMinimumLeafSize());
		newTree->SetLeftNode(NewLeft); newTree->SetRightNode(NewRight);

		// Set splitting variables and splitting rules.
		// Need to compute the potential.
		int newVar = ran.Discrete(prob);
		double newTau = dens[newVar]->Sample(ran);
		if (newVar<0)
			Rcout<<"Something wrong in ProposalPruneGrow: newVar is faulty!";
		//int oldVar = newTree->GetSplitVariable();
		//double oldTau = newTree->GetSplitLevel();

		newTree->SetSplitVariable(newVar);
		newTree->SetSplitLevel(newTau);

		potential = 0.0;
		// Here the probability to grow a tree is 1, while the reverse has a prob. of 0.5 to prune.
		potential -= -log((double) prob[newVar]);
		potential -= dens[newVar]->Potential(newTau,ran);
		potential += -log(0.5);
		// Update Subject List

		if (!newTree->UpdateSubjectList())
		{
			Delta *delta=NULL;
			newTree->RemoveSubTree();
			delete newTree;
			newTree=NULL;
			return delta;
		}
		else
		{
			Delta *delta = new DeltaGrow(&tree,newVar,newTau);
			newTree->RemoveSubTree();
			delete newTree;
			newTree=NULL;
			return delta;
		};
	};



	//
	// Now we know the root is interior. Find all interior nodes
	//

	double coin=ran.Unif01();

	if (coin<=0.5)
		// In the case where grow move is refrained, we need to change this value to larger one so that grow move will be encouraged.
		// || tree.GetSize()==10) // Just for test CAUTION!!!! Don't leave this code here!!!!
	{
		// To grow the tree
		//
		// draw one of the leaf nodes uniformly
		//
		std::vector<const Node *> node;
	  std::vector<const Node *> leaves;

		node.push_back(&tree);
		int i;
		for (i = 0; i < node.size(); i++)
		{
			if (node[i]->GetRightNode()==NULL)
				leaves.push_back(node[i]);
			else
			{
				node.push_back(node[i]->GetRightNode());
				node.push_back(node[i]->GetLeftNode());
			};
		};



		std::vector<double> pr = std::vector<double>(leaves.size(),1.0/leaves.size());
		int nr = ran.Discrete(pr);  // i.e. propose new rule for this node

		Node *NewNode = leaves[nr]->Copy();
		Node *NewLeft=new Node;
		NewLeft->SetObservation(leaves[nr]->GetObservation());
		NewLeft->SetSplitVariable(-1);
		NewLeft->SetMinimumLeafSize(leaves[nr]->GetMinimumLeafSize());
		Node *NewRight=new Node;
		NewRight->SetSplitVariable(-1);
		NewRight->SetObservation(leaves[nr]->GetObservation());
		NewRight->SetMinimumLeafSize(leaves[nr]->GetMinimumLeafSize());

		NewNode->SetLeftNode(NewLeft);
		NewNode->SetRightNode(NewRight);

		// Need to set the splitting variables and threshold.
		node.clear();
		node.push_back(&tree);
		std::vector<const Node *> Possible;

		for (i = 0; i < node.size(); i++)
		{
			if (node[i]->GetRightNode()!=NULL)
			{
				if (node[i]->GetRightNode()->GetRightNode()==NULL && node[i]->GetLeftNode()->GetRightNode()==NULL &&
					node[i]->GetRightNode() != leaves[nr] && node[i]->GetLeftNode() != leaves[nr])
					Possible.push_back(node[i]); // These are the possible node to be pruned
				node.push_back(node[i]->GetRightNode());
				node.push_back(node[i]->GetLeftNode());
			};
		};
		Possible.push_back(NULL);  // the new interior node made when groing the tree.

		if (Possible.size()==0)
		{
			Rcout<<"No node can be pruned. Something is wrong."<<endl;
		};

		int newVar = ran.Discrete(prob);
		double newTau = dens[newVar]->Sample(ran);
		NewNode->SetSplitVariable(newVar);
		NewNode->SetSplitLevel(newTau);


		if (newVar<0)
		  Rcout<<"Something wrong in ProposalPruneGrow: newVar is faulty!";

		potential = 0.0;
		potential -= -log((double) prob[newVar]);
		potential -= dens[newVar]->Potential(newTau,ran);
		potential -= -log((double) 1.0/leaves.size());
		potential += -log((double) 1.0/Possible.size());

		//
		// Then we have to update the subject list for all
		// descendent of the chosen node. Let a function
		// in the Node do this.
		//

		if (!NewNode->UpdateSubjectList())
		{
			Delta *delta=NULL;
			NewNode->RemoveSubTree();
			delete NewNode;
			NewNode=NULL;

			return delta;
		}
		else
		{
			Delta *delta = new DeltaGrow(leaves[nr],newVar,newTau);
			NewNode->RemoveSubTree();
			delete NewNode;
			NewNode = NULL;

			return delta;
		};
	}
	else
	{
		// To prune a tree
		std::vector<const Node *> node;
	  std::vector<const Node *> Possible;
		node.push_back(&tree);
		int i;
		for (i = 0; i < node.size(); i++)
		{
			if (node[i]->GetRightNode()!=NULL)
			{
				if (node[i]->GetRightNode()->GetRightNode()==NULL && node[i]->GetLeftNode()->GetRightNode()==NULL)
					Possible.push_back(node[i]); // These are the possible node to be pruned
				node.push_back(node[i]->GetRightNode());
				node.push_back(node[i]->GetLeftNode());
			};
		};

		if (Possible.size()==0)
		{
			Rcout<<"No node can be pruned. Something is wrong."<<endl;
		};

		// Uniformly pick a possible node to prune
		std::vector<double> pr = std::vector<double>(Possible.size(),1.0/Possible.size());
		int nr = ran.Discrete(pr);
		int oldVar=Possible[nr]->GetSplitVariable();
		double oldTau=Possible[nr]->GetSplitLevel();

		node.clear();
		std::vector<const Node *> leaves;
		node.push_back(&tree);

		for (i = 0; i < node.size(); i++)
		{
			if (node[i]->GetRightNode()==NULL)
				leaves.push_back(node[i]);
			else
			{
				node.push_back(node[i]->GetRightNode());
				node.push_back(node[i]->GetLeftNode());
			};
		};

		potential = 0;
		potential += -log((double) prob[oldVar]);
		potential += dens[oldVar]->Potential(oldTau,ran);
		potential += -log((double) 1.0/(leaves.size()-1));
		potential -= -log((double) 1.0/Possible.size());

		// Hakon suggest, in the case the tree is about to be pruned to a root, subtract -log(0.5);

		if (Possible.size()==1 && leaves.size()==2) potential -= -log(0.5);

		// About newvar and newtau here;

		Delta *delta = new DeltaPrune(Possible[nr]);
		return delta;

	}; // end of if-else


	//  NodeTree *newTree = tree.CopyTree();

	Delta *delta = NULL;
	potential = 0.0;

	return delta;
};
