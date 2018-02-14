#include <cmath>

#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "NodeTree.h"
#include "Delta.h"
#include "DP.h"
#include "ModelTreeStructurePinBall.h"


ModelTreeStructurePinBall::ModelTreeStructurePinBall(double paramNumberOfLeaves,
													 double paramSplit) :
ModelTreeStructure()
{
	lambda = paramNumberOfLeaves;
	prob = paramSplit;
	MinimumSize = 1;

	return;
};

ModelTreeStructurePinBall::ModelTreeStructurePinBall(double paramNumberOfLeaves,
													 double paramSplit, int MinimumTreeSize) :
ModelTreeStructure()
{
	lambda = paramNumberOfLeaves;
	prob = paramSplit;
	MinimumSize = MinimumTreeSize;

	return;
};


ModelTreeStructurePinBall::~ModelTreeStructurePinBall(void)
{
	return;
};


double ModelTreeStructurePinBall::PotentialDifference(const NodeTree &oldTree,
													  const Delta &newTree,
													  Random &ran) const
{
	double pot = DP::PotentialDifference(oldTree,newTree,*this,ran);

	return pot;
};



NodeTree *ModelTreeStructurePinBall::Simulate(Random &ran, Observation *obs) const
{
	int nLeaf = 1 + ran.Poisson(lambda);
	NodeTree *root = new NodeTree(obs);

	std::vector<Node *> node;
	std::vector<int> size;
	node.push_back(root);
	size.push_back(nLeaf);

	int i;
	for (i = 0; i < node.size(); i++)
	{
		if (size[i] > 1)
		{
			int nRight;
			if (ran.Unif01() <= 0.5)
				nRight = 1 + ran.Binomial(size[i]-2,prob);
			else
				nRight = 1 + ran.Binomial(size[i]-2,1.0-prob);
			int nLeft = size[i] - nRight;

			Node *right = new Node;
			Node *left = new Node;
			right->SetObservation(obs);
			right->SetParentNode(node[i]);
			left->SetObservation(obs);
			left->SetParentNode(node[i]);
			node[i]->SetRightNode(right);
			node[i]->SetLeftNode(left);

			node.push_back(right);
			node.push_back(left);
			size.push_back(nRight);
			size.push_back(nLeft);
		}
	}

	return root;
};


NodeTree *ModelTreeStructurePinBall::Simulate(Random &ran, Observation *obs, int LeafSize) const
{
	int nLeaf = MinimumSize + ran.Poisson(lambda);
	//vector <double> nLeafProbs=GetnLeafProbs(floor(1.0*obs->GetN()/LeafSize), ran);
	//vector <double> nLeafProbs=GetnLeafProbs(6, ran);
	//int nLeaf=ran.Discrete(nLeafProbs)+1;
	//int nLeaf=1; //Starting from root;
	NodeTree *root = new NodeTree(obs);
	int i;
	Rcout<<"Number of leaves: "<<nLeaf<<endl;
	std::vector<Node *> node;
	std::vector<int> size;
	node.push_back(root);
	size.push_back(nLeaf);
	//int i;
	for (i = 0; i < node.size(); i++)
	{
		if (size[i] > 1)
		{
			int nRight;
			if (ran.Unif01() <= 0.5)
				nRight = 1 + ran.Binomial(size[i]-2,prob);
			else
				nRight = 1 + ran.Binomial(size[i]-2,1.0-prob);
			int nLeft = size[i] - nRight;

			Node *right = new Node;
			Node *left = new Node;
			right->SetObservation(obs);
			right->SetParentNode(node[i]);
			right->SetMinimumLeafSize(LeafSize);
			left->SetObservation(obs);
			left->SetParentNode(node[i]);
			left->SetMinimumLeafSize(LeafSize);
			node[i]->SetRightNode(right);
			node[i]->SetLeftNode(left);

			node.push_back(right);
			node.push_back(left);
			size.push_back(nRight);
			size.push_back(nLeft);
		};
	};

	return root;
};



Model *ModelTreeStructurePinBall::Copy(void) const
{
	Model *model = new ModelTreeStructurePinBall(lambda,prob, MinimumSize);

	return model;
};



double ModelTreeStructurePinBall::Potential(const NodeTree &tree,Random &ran) const
{
  std::vector<double> nNode;
  std::vector<double> left;
  std::vector<double> right;
  std::vector<double> parent;
  std::vector<const Node *> node;

	node.push_back(&tree);
	nNode.push_back(0);
	left.push_back(-1);
	right.push_back(-1);
	parent.push_back(-1);

	int i;
	for (i = 0; i < node.size(); i++)
	{
		const Node *ll = node[i]->GetLeftNode();
		const Node *rr = node[i]->GetRightNode();

		if (ll != NULL)
		{
			node.push_back(ll);
			nNode.push_back(-1);
			left.push_back(-1);
			right.push_back(-1);
			parent.push_back(i);

			left[i] = left.size() - 1;
		}

		if (rr != NULL)
		{
			node.push_back(rr);
			nNode.push_back(-1);
			left.push_back(-1);
			right.push_back(-1);
			parent.push_back(i);

			right[i] = right.size() - 1;
		}
	}

	for (i = node.size() - 1; i >= 0; i--)
	{
		if (left[i] == -1 && right[i] == -1)
			nNode[i] = 1;
		else
			nNode[i] = nNode[left[i]] + nNode[right[i]];
	}

	double pot=0.0;
	int nLeaf=nNode[0];
	//  pot=ran.PotentialTruncatedPoisson(lambda, nLeaf, floor(1.0*tree.GetObservation()->GetN()/tree.GetMinimumLeafSize()));

	if (nLeaf<MinimumSize)
	  Rcout<<"Less than the minimum leaf number!"<<endl;
	pot = ran.PotentialPoisson(lambda,nLeaf-MinimumSize);

	for (i = 0; i < node.size(); i++)
	{
		if (nNode[i] > 1)
		{
			int nRight = nNode[right[i]];
			double pot1 = ran.PotentialBinomial(nNode[i]-2,prob,nRight-1);
			double pot2 = ran.PotentialBinomial(nNode[i]-2,1.0-prob,nRight-1);

			pot += - log(0.5*exp(- pot1) + 0.5*exp(- pot2)); //Miss 1/2?
		}
	}

	return pot;
};

double ModelTreeStructurePinBall::PotentialDifferencePrune(const NodeTree &oldTree,const Node *newLeaf,Random &ran) const
{
  std::vector<double> nNode;
  std::vector<double> left;
  std::vector<double> right;
  std::vector<double> parent;
  std::vector<const Node *> node;

	int nrNewLeaf = -1;

	node.push_back(&oldTree);
	nNode.push_back(0);
	left.push_back(-1);
	right.push_back(-1);
	parent.push_back(-1);

	int i;
	for (i = 0; i < node.size(); i++)
	{
		if (node[i] == newLeaf) nrNewLeaf = i;

		const Node *ll = node[i]->GetLeftNode();
		const Node *rr = node[i]->GetRightNode();

		if (ll != NULL)
		{
			node.push_back(ll);
			nNode.push_back(-1);
			left.push_back(-1);
			right.push_back(-1);
			parent.push_back(i);

			left[i] = left.size() - 1;
		}

		if (rr != NULL)
		{
			node.push_back(rr);
			nNode.push_back(-1);
			left.push_back(-1);
			right.push_back(-1);
			parent.push_back(i);

			right[i] = right.size() - 1;
		}
	}

	for (i = node.size() - 1; i >= 0; i--)
	{
		if (left[i] == -1 && right[i] == -1)
			nNode[i] = 1;
		else
			nNode[i] = nNode[left[i]] + nNode[right[i]];
	}


	double pot=0.0;
	int nLeaf=nNode[0];
	//  pot = - ran.PotentialTruncatedPoisson(lambda, nLeaf, floor(1.0*oldTree.GetObservation()->GetN()/oldTree.GetMinimumLeafSize()));
	if (nLeaf<MinimumSize)
	  Rcout<<"Less than the minimum leaf number!"<<endl;

	pot = - ran.PotentialPoisson(lambda,nLeaf-MinimumSize);

	for (i = 0; i < node.size(); i++)
	{
		if (nNode[i] > 1)
		{
			int nRight = nNode[right[i]];
			double pot1 = ran.PotentialBinomial(nNode[i]-2,prob,nRight-1);
			double pot2 = ran.PotentialBinomial(nNode[i]-2,1.0-prob,nRight-1);

			pot -= - log(0.5*exp(- pot1) + 0.5*exp(- pot2)); //Miss 1/2?
		}
	}


	int pp;
	nNode[nrNewLeaf] = 1;
	for (pp = parent[nrNewLeaf]; pp >= 0; pp = parent[pp])
		nNode[pp] -= 1;

	nLeaf = nNode[0];
	//  pot += ran.PotentialTruncatedPoisson(lambda, nLeaf, floor(1.0*oldTree.GetObservation()->GetN()/oldTree.GetMinimumLeafSize()));

	if (nLeaf<MinimumSize)
	  Rcout<<"Less than the minimum leaf number!"<<endl;
	pot += ran.PotentialPoisson(lambda,nLeaf-MinimumSize);

	for (i = 0; i < node.size(); i++)
	{
		if (nNode[i] > 1)
		{
			int nRight = nNode[right[i]];
			double pot1 = ran.PotentialBinomial(nNode[i]-2,prob,nRight-1);
			double pot2 = ran.PotentialBinomial(nNode[i]-2,1.0-prob,nRight-1);

			pot += - log(0.5*exp(- pot1) + 0.5*exp(- pot2)); //Miss 1/2?
		}
	}

	return pot;
};






double ModelTreeStructurePinBall::PotentialDifferenceGrow(const NodeTree &oldTree,const Node *oldLeaf,Random &ran) const
{
  std::vector<double> nNode;
  std::vector<double> left;
  std::vector<double> right;
  std::vector<double> parent;
  std::vector<const Node *> node;

	int nrOldLeaf = -1;

	node.push_back(&oldTree);
	nNode.push_back(0);
	left.push_back(-1);
	right.push_back(-1);
	parent.push_back(-1);

	int i;
	for (i = 0; i < node.size(); i++)
	{
		if (node[i] == oldLeaf) nrOldLeaf = i;

		const Node *ll = node[i]->GetLeftNode();
		const Node *rr = node[i]->GetRightNode();

		if (ll != NULL)
		{
			node.push_back(ll);
			nNode.push_back(-1);
			left.push_back(-1);
			right.push_back(-1);
			parent.push_back(i);

			left[i] = left.size() - 1;
		}

		if (rr != NULL)
		{
			node.push_back(rr);
			nNode.push_back(-1);
			left.push_back(-1);
			right.push_back(-1);
			parent.push_back(i);

			right[i] = right.size() - 1;
		}
	}

	for (i = node.size() - 1; i >= 0; i--)
	{
		if (left[i] == -1 && right[i] == -1)
			nNode[i] = 1;
		else
			nNode[i] = nNode[left[i]] + nNode[right[i]];
	}

	double pot=0.0;
	int nLeaf=nNode[0];
	//  pot = - ran.PotentialTruncatedPoisson(lambda, nLeaf, floor(1.0*oldTree.GetObservation()->GetN()/oldTree.GetMinimumLeafSize()));

	if (nLeaf<MinimumSize)
	  Rcout<<"Less than the minimum leaf number!"<<endl;
	pot = - ran.PotentialPoisson(lambda,nLeaf-MinimumSize);

	for (i = 0; i < node.size(); i++)
	{
		if (nNode[i] > 1)
		{
			int nRight = nNode[right[i]];
			double pot1 = ran.PotentialBinomial(nNode[i]-2,prob,nRight-1);
			double pot2 = ran.PotentialBinomial(nNode[i]-2,1.0-prob,nRight-1);

			pot -= - log(0.5*exp(- pot1) + 0.5*exp(- pot2)); //Miss 1/2?
		}
	}


	nNode.push_back(1);
	left.push_back(-1);
	right.push_back(-1);
	parent.push_back(nrOldLeaf);
	left[nrOldLeaf] = left.size() - 1;

	nNode.push_back(1);
	left.push_back(-1);
	right.push_back(-1);
	parent.push_back(nrOldLeaf);
	right[nrOldLeaf] = right.size() - 1;

	int pp;
	nNode[nrOldLeaf] = 2;
	for (pp = parent[nrOldLeaf]; pp >= 0; pp = parent[pp])
		nNode[pp] += 1;

	nLeaf = nNode[0];
	//  pot += ran.PotentialTruncatedPoisson(lambda, nLeaf, floor(1.0*oldTree.GetObservation()->GetN()/oldTree.GetMinimumLeafSize()));

	if (nLeaf<MinimumSize)
	  Rcout<<"Less than the minimum leaf number!"<<endl;
	pot += ran.PotentialPoisson(lambda,nLeaf-MinimumSize);

	for (i = 0; i < node.size(); i++)
	{
		if (nNode[i] > 1)
		{
			int nRight = nNode[right[i]];
			double pot1 = ran.PotentialBinomial(nNode[i]-2,prob,nRight-1);
			double pot2 = ran.PotentialBinomial(nNode[i]-2,1.0-prob,nRight-1);

			pot += - log(0.5*exp(- pot1) + 0.5*exp(- pot2)); //Miss 1/2?
		}
	}

	return pot;
};


std::vector <double> ModelTreeStructurePinBall::GetnLeafProbs(int length, Random &ran) const
{
  std::vector <double> Probs(length);
	int i;
	double sum=0;
	for (i=0;i<length;i++)
		Probs[i]=exp(-ran.PotentialPoisson(lambda, i+1));
	// Normalize Probs;
	for (i=0;i<length;i++)
		sum+=Probs[i];
	for (i=0;i<length;i++)
		sum=Probs[i]/sum;
	return Probs;
};


