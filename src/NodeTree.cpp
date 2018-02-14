#include <sstream>
#include "Random.h"
#include "Model.h"
#include "Observation.h"
#include "Density.h"
#include "NodeTree.h"
#include <numeric>

NodeTree::NodeTree()
{
	parent=NULL;
	this->nLeaf=0;
	this->nLevel=0;
	this->splitVariable=-1;
	left=NULL;
	right=NULL;
	return;
};

NodeTree::NodeTree(const Observation *obs)
{
	this->SetObservation(obs);
	int n=this->GetObservation()->GetN();
	this->SubjectList.resize(n);
	for (int i=0;i<n;i++)
		SubjectList[i]=i; //Initially, SubjectList has all the observations in the dataset.*/
	parent=NULL;
	left=NULL;
	right=NULL;
	splitVariable=-1;
	for (int i = 0; i < obs->GetP(); i++)
	  Splitting_var.push_back(i);
	this->nLeaf=0;
	this->nLevel=0;
	return;
};

NodeTree::~NodeTree()
{
	if (left!=NULL)
	{
		left->RemoveSubTree();
		delete left;
		left=NULL;
	};
	if (right!=NULL)
	{
		right->RemoveSubTree();
		delete right;
		right=NULL;
	};
	return;
};

NodeTree::NodeTree(Node &root)
{
	parent=NULL;
	this->nLeaf=0;
	this->nLevel=0;
	Node* NewNode=root.CopySubTree();
	this->SetLeftNode(NewNode->GetLeftNode());
	this->SetRightNode(NewNode->GetRightNode());
	this->SetObservation(NewNode->GetObservation());
	this->SetSplitLevel(NewNode->GetSplitLevel());
	this->SetSplitVariable(NewNode->GetSplitVariable());
	this->SetSubjectList(NewNode->GetSubjectList());
	this->SetSplittingVar(NewNode->GetSplittingVar());
	this->SetMinimumLeafSize(NewNode->GetMinimumLeafSize());
	return;
}

void NodeTree::SetParentNode(Node *node)
{
	Rcout<<"The parent for a tree can only be NULL"<<endl;
};

NodeTree* NodeTree::CopyTree(void) const
{
	NodeTree* NewTree=new NodeTree;
	NewTree->SetObservation(this->GetObservation());
	NewTree->SetSplitLevel(this->GetSplitLevel());
	NewTree->SetSplitVariable(this->GetSplitVariable());
	NewTree->SetSplittingVar(this->GetSplittingVar());
	NewTree->SetSubjectList(this->GetSubjectList());
	NewTree->SetMinimumLeafSize(this->GetMinimumLeafSize());
	if (this->GetLeftNode()!=NULL || this->GetRightNode()!=NULL)
	{
		Node* pLeft=this->GetLeftNode()->CopySubTree();
		Node* pRight=this->GetRightNode()->CopySubTree();
		NewTree->SetLeftNode(pLeft);
		NewTree->SetRightNode(pRight);
	}

	return NewTree;
};

double NodeTree::Pot(const Model &modelStructure, const Model &modelVar,  const Model &modelLikelihood, Random &ran)
{
	double pot=0;
	pot =  modelStructure.Potential(*this,ran);
	pot += modelVar.Potential(*this, ran);
	pot += modelLikelihood.Potential(*this, ran);
	potential=pot;
	return pot;
};

double NodeTree::LogLikelihood(const Model &modelLikelihood, Random &ran)
{
	double pot=0;
	pot += modelLikelihood.Potential(*this, ran);
	return -pot;
};

int NodeTree::GetSize() const
{
	NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;
	int Size;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};
	Size=leaves.size();
	delete Tree;

	return Size;
};

std::vector <Node *> NodeTree::GetLeaves()
{
	NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};
	return leaves;
};

List NodeTree::GetVal()
{
  double value = 0;
  std::vector <int> label;
  
  NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
  node.push_back(Tree);
  std::vector <Node *> leaves;
  int i;
  
  for (i=0;i<node.size();i++)
  {
    if (node[i]->GetRightNode()==NULL)
    {
      leaves.push_back(node[i]);
    }
    else
    {
      node.push_back(node[i]->GetLeftNode());
      node.push_back(node[i]->GetRightNode());
    };
  };
  
  for (i=0;i<leaves.size();i++) //loop through leaves
  {
    std::vector <int> SubjectList=leaves[i]->GetSubjectList();
    
    double max = 0;
    int opt = 0;
    for (int k = 0; k < (obs->GetNumCats()+1); k++) //loop through treatments
    {
      double temp = 0;
      for (int j=0; j<SubjectList.size(); j++){
        temp = temp + obs->GetV(SubjectList[j], k);
      };
      if (temp > max) {
        opt = k;
        max = temp;
      };
    };
    label.push_back(opt);
    value = value + max;
  };
  
  delete Tree;
  leaves.clear();
  node.clear();
  
  return List::create(Named("label")=label,
                      Named("value")=value);
};

std::vector <std::vector <double> > NodeTree::GetLeavesObservations()
{
  std::vector <std::vector <double> > Result;
	//vector <Node *> leaves=GetLeaves();
	//int i;

	NodeTree *Tree=this->CopyTree();
	std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};

	for (i=0;i<leaves.size();i++)
	{
	  std::vector <double> temp;
	  std::vector <int> SubjectList=leaves[i]->GetSubjectList();
		for (int j=0;j<SubjectList.size();j++)
			temp.push_back(obs->GetY(SubjectList[j]));
		Result.push_back( temp );
	};
	delete Tree;
	leaves.clear();
	node.clear();


	return Result;
};

std::vector <std::vector <int> > NodeTree::GetLeavesSubjects()
{
  std::vector <std::vector <int> > Result;
	//vector <Node *> leaves=GetLeaves();
	//int i;

	NodeTree *Tree=this->CopyTree();
	std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};

	for (i=0;i<leaves.size();i++)
	{
	  std::vector <int> SubjectList=leaves[i]->GetSubjectList();
		Result.push_back(SubjectList);
	};
	delete Tree;
	leaves.clear();
	node.clear();

	return Result;
};

std::vector <int> NodeTree::GetSplittingVariableSet() const
{
  std::vector <int> Result;
	//vector <Node *> leaves=GetLeaves();
	//int i;

	NodeTree *Tree=this->CopyTree();
	std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{
			Result.push_back(node[i]->GetSplitVariable());
			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};

	delete Tree;
	leaves.clear();
	node.clear();

	return Result;
};

bool NodeTree::Misclassfication(int &LeavesSize, int &Error, Density *like) const
{

	NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;
	int k = Tree->GetObservation()->GetK();

	LeavesSize=0; Error=0;


	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL) //leaf
		{
			leaves.push_back(node[i]);

		  std::vector<int> subject(node[i]->GetSubjectList());
			const Observation *obs = node[i]->GetObservation();
			std::vector <int> s1;

			int l = subject.size();

			int s,j;

			NumericVector leaf;

			for (s = 0; s < l; s++)
			{
			  int nr = subject[s];
			  leaf.push_back(obs->GetY(nr));
			  for (j = 0; j < k; j++){
			    int count = 0;
			    if (obs->GetY(nr)==j) count++;
			    s1.push_back(count);
			  };
			};

			std::vector <double> pMean=like->PosteriorMean(leaf);

			pMean.push_back(1-std::accumulate(pMean.begin(), pMean.end(), 0.0));
			int pmax=std::distance(pMean.begin(), std::max_element(pMean.begin(), pMean.end()));
			if (pmax==k)
			  s1.erase(s1.begin());
			else
			  s1.erase(s1.begin()+pmax-1);

			Error+=std::accumulate(s1.begin(), s1.end(), 0);
		}
		else
		{
			node.push_back(node[i]->GetRightNode());
			node.push_back(node[i]->GetLeftNode());
		};
	};
	LeavesSize=leaves.size();
	delete Tree;
	Tree=NULL;
	return true;
};



