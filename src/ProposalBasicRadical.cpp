#include "Random.h"
#include "Observation.h"
#include "Model.h"
#include "Density.h"
#include "NodeTree.h"
#include "Delta.h"
#include "DeltaBasic.h"
#include "ProposalBasicRadical.h"
#include "ChangeLeaves.h"

ProposalBasicRadical::ProposalBasicRadical(const ChangeLeaves &changeLeaves)
{
	this->changeLeaves = changeLeaves.copy();

	return;
};


ProposalBasicRadical::~ProposalBasicRadical(void)
{
	delete changeLeaves;

	return;
};


std::vector <std::vector <int> > ProposalBasicRadical::GetLeavesConfiguration(const NodeTree &tree, std::vector <int> &CurrentSplittingVariables) const
{
	//	NodeTree *Tree=(&tree)->CopySubTree();
	NodeTree *newTree = tree.CopyTree();
  std::vector <Node *> node;
  std::vector <Node *> leaves;
  std::vector <std::vector <int> > Configuration;
	int i;

	node.push_back(newTree);
	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{
			node.push_back(node[i]->GetRightNode());
			node.push_back(node[i]->GetLeftNode());
		};
	};
	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetSplitVariable()!=-1)
			CurrentSplittingVariables.push_back(node[i]->GetSplitVariable());
	};

	for (i=0;i<leaves.size();i++)
	{
	  std::vector<int> list(leaves[i]->GetSubjectList());
		Configuration.push_back(list);
	};

	//ShowLeaves(Configuration);
	newTree->RemoveSubTree();
	delete newTree;
	newTree=NULL;
	return Configuration;
};


Delta * ProposalBasicRadical::ProposeChange(const NodeTree &tree, double &potential, Random &ran) const
{

  std::vector <int> CurrentSplittingVariables;
  std::vector <std::vector <int> > Leaves=GetLeavesConfiguration(tree, CurrentSplittingVariables);

	double pot=0;
	pot -= changeLeaves->ProposeChange(Leaves,tree.GetObservation(),tree.GetMinimumLeafSize(),ran);

	std::vector <int> NewSplittingVariables;
	bool STOP=false;
	const Observation *obs=tree.GetObservation();

	std::vector <Node *> node;
	NodeTree* NewTree=new NodeTree(obs);
	std::vector <std::vector <std::vector <int> > > LeavesList;

	NewTree->SetMinimumLeafSize(tree.GetMinimumLeafSize());

	//ShowLeaves(Leaves);

	if (Leaves.empty())    // no interior nodes exist!
	{
		Rcout<<"No interior node in restructure. Terminated.\n";
	};

	if (Leaves.size()==1)
	{
		NewTree->RemoveSubTree();
		delete NewTree;
		NewTree=NULL;
		Delta *delta = NULL;
		return delta;
	};

	//ExtendSplittingVariables(CurrentSplittingVariables, NewSplittingVariables, pot, true);
	for (int i=0;i<obs->GetP();i++)
		NewSplittingVariables.push_back(i);
	// The new splitting variable will be all the possible predictors.

	//ShowSplittingVariables(NewSplittingVariables);
	node.push_back(NewTree);
	LeavesList.push_back(Leaves);

	for(int i=0;i<node.size();i++)
	{
	  std::vector <int> SubjectList=node[i]->GetSubjectList();
	  std::vector <int> SplittingVariables;
	  std::vector <double> SplittingThreshold;
	  std::vector <double> SplittingThresholdProb;

		Leaves=LeavesList[i];

		if (Leaves.size()>1)
		{
			for (int j=0;j<NewSplittingVariables.size();j++) //figure 3.8: finding possible intervals
			{
				int PossibleSplittingVariables=NewSplittingVariables[j];
			  std::vector <double> max(Leaves.size());//max and min are vectors of maxs and mins of j-th var in each leaf
			  std::vector <double> min(Leaves.size());

				//ShowLeaves(Leaves);
				for (int k=0;k<max.size();k++)
				{
					max[k]=obs->GetX(Leaves[k][0], PossibleSplittingVariables);
					min[k]=obs->GetX(Leaves[k][0], PossibleSplittingVariables);

					for (int l=1;l<Leaves[k].size();l++)
					{
						if (obs->GetX(Leaves[k][l], PossibleSplittingVariables) >max[k])
							max[k]=obs->GetX(Leaves[k][l], PossibleSplittingVariables);
						if (obs->GetX(Leaves[k][l], PossibleSplittingVariables)<min[k])
							min[k]=obs->GetX(Leaves[k][l], PossibleSplittingVariables);
					};
				};

				Sort(max, min); // Sort max ascendently

				for (int k=0;k<max.size()-1;k++)
				{
					bool Possible=true;

					double CloseMin=max[max.size()-1];
					for (int l=k+1;l<max.size() && Possible;l++)
					{
						if (min[l]<=max[k]) //overlapping, no way to set splitting here
						{
							Possible=false;
						};
						if (min[l]>max[k] && min[l]<CloseMin)
							CloseMin=min[l];
					};
					if (Possible)
					{
						double threshold=ran.Unif01()*(CloseMin-max[k])+max[k];
						// Proposed a uniform distribution as new threshold so far;
						// Could use diffrent distribution here.

						SplittingVariables.push_back(PossibleSplittingVariables);
						SplittingThreshold.push_back(threshold);
						SplittingThresholdProb.push_back(CloseMin-max[k]);
					};
				};
			};

			if (SplittingVariables.empty())
			{
				NewTree->RemoveSubTree();
				delete NewTree;
				NewTree=NULL;
				Delta *delta=NULL;
				return delta;
			};

			// Choose one of the possible splits uniformly
			std::vector <double> Probs(SplittingVariables.size(), 1.0/double(SplittingVariables.size()));
			if (SplittingVariables.size()!=SplittingThreshold.size()){
			  Rcout<<"Something wrong with the splitting variables!"<<endl;
			};

			int Pick=ran.Discrete(Probs);
			int NewVariable=SplittingVariables[Pick];
			double NewThreshold=SplittingThreshold[Pick];

			pot-= -log(1.0/SplittingVariables.size());
			pot-= -log(1.0/SplittingThresholdProb[Pick]);

			node[i]->SetSplitVariable(NewVariable);
			node[i]->SetSplitLevel(NewThreshold);

			Node *Left=new Node;
			Node *Right=new Node;

			std::vector <std::vector <int> > NewLeftLeaves, NewRightLeaves;
			std::vector <int> NewLeftSubjectList, NewRightSubjectList;

			for (int j=0;j<Leaves.size();j++)
			{
			  if (obs->GetX(Leaves[j][0], NewVariable) <= NewThreshold)
			    NewLeftLeaves.push_back(Leaves[j]);
			  if (obs->GetX(Leaves[j][0], NewVariable) > NewThreshold)
			    NewRightLeaves.push_back(Leaves[j]);
			};

			for (int j=0;j<SubjectList.size();j++)
			{
			  if (obs->GetX(SubjectList[j], NewVariable) <= NewThreshold)
			    NewLeftSubjectList.push_back(SubjectList[j]);
			  if (obs->GetX(SubjectList[j], NewVariable) > NewThreshold)
			    NewRightSubjectList.push_back(SubjectList[j]);
			};

			Right->SetObservation(obs);
			Left->SetObservation(obs);
			Right->SetParentNode(node[i]);
			Left->SetParentNode(node[i]);
			Right->SetSubjectList(NewRightSubjectList);
			Left->SetSubjectList(NewLeftSubjectList);
			Right->SetMinimumLeafSize(node[i]->GetMinimumLeafSize());
			Left->SetMinimumLeafSize(node[i]->GetMinimumLeafSize());
			node[i]->SetLeftNode(Left);
			node[i]->SetRightNode(Right);
			node.push_back(Left);
			node.push_back(Right);

			LeavesList.push_back(NewLeftLeaves);
			LeavesList.push_back(NewRightLeaves);
		};
	};

	//	Compute P(x|r,y)*P(r|y)
	//	*****************************************************

	//CurrentSplittingVariables.clear();
	Leaves.clear();
	std::vector <int> Test;
	Leaves=GetLeavesConfiguration(*NewTree, Test);

	pot += changeLeaves->ProposeReverse(Leaves,NewTree->GetObservation(),NewTree->GetMinimumLeafSize(),ran);

	STOP=false;

	node.clear();

	if (Leaves.empty())    // no interior nodes exist!
	{
		Rcout<<"no interior nodes exist! Terminated.\n";
	};

	node.push_back(NewTree);
	LeavesList.push_back(Leaves);

	for(int i=0;i<node.size();i++)
	{
	  std::vector <int> SubjectList=node[i]->GetSubjectList();
	  std::vector <int> SplittingVariables;
	  std::vector <double> SplittingThreshold;
	  std::vector <double> SplittingThresholdProb;

		Leaves=LeavesList[i];
		//ShowSubjectList(SubjectList);

		if (Leaves.size()>1)
		{
		  for (int j=0;j<NewSplittingVariables.size();j++)
			{
				int PossibleSplittingVariables=NewSplittingVariables[j];
			  std::vector <double> max(Leaves.size());
			  std::vector <double> min(Leaves.size());

				//ShowLeaves(Leaves);
				for (int k=0;k<max.size();k++)
				{
					max[k]=obs->GetX(Leaves[k][0], PossibleSplittingVariables);
					min[k]=obs->GetX(Leaves[k][0], PossibleSplittingVariables);

					for (int l=1;l<Leaves[k].size();l++)
					{
						if (obs->GetX(Leaves[k][l], PossibleSplittingVariables) >max[k])
							max[k]=obs->GetX(Leaves[k][l], PossibleSplittingVariables);
						if (obs->GetX(Leaves[k][l], PossibleSplittingVariables)<min[k])
							min[k]=obs->GetX(Leaves[k][l], PossibleSplittingVariables);
					};
				};

				Sort(max, min); // Sort max ascendently

				for (int k=0;k<max.size()-1;k++)
				{
					bool Possible=true;

					double CloseMin=max[max.size()-1];
					for (int l=k+1;l<max.size() && Possible;l++)
					{
						if (min[l]<=max[k])
						{
							Possible=false;
						};
						if (min[l]>max[k] && min[l]<CloseMin)
							CloseMin=min[l];
					};
					if (Possible)
					{
						double threshold=ran.Unif01()*(CloseMin-max[k])+max[k];
						// Proposed a uniform distribution as new threshold so far;
						// Could use diffrent distribution here.
						if (threshold>=1)
						{
							Rcout<<"Sth. wrong with threshold...\n";
						};
						SplittingVariables.push_back(PossibleSplittingVariables);
						SplittingThreshold.push_back(threshold);
						SplittingThresholdProb.push_back(CloseMin-max[k]);
					};
				};
			};

			if (SplittingVariables.empty())
			{
				Rcout<<"Something Wrong. Terminated";
				NewTree->RemoveSubTree();
				delete NewTree;
				NewTree=NULL;
				Delta *delta=NULL;
				return delta;
			};

			// Choose one of the possible splits uniformly
			std::vector <double> Probs(SplittingVariables.size(), 1.0/double(SplittingVariables.size()));
			if (SplittingVariables.size()!=SplittingThreshold.size()){
			  Rcout<<"Something wrong with the splitting variables!"<<endl;
			}

			int TrueSplittingVariable=node[i]->GetSplitVariable();

			int Pick=0; //ran.Discrete(Probs);
			while (Pick<SplittingVariables.size())
			{
				if (SplittingVariables[Pick]==TrueSplittingVariable)
					break;
				Pick++;
			};

			int NewVariable=SplittingVariables[Pick];
			double NewThreshold=node[i]->GetSplitLevel();

			pot+= -log(1.0/SplittingVariables.size());
			pot+= -log(1.0/SplittingThresholdProb[Pick]);

			//node[i]->SetSplitVariable(NewVariable);
			//node[i]->SetSplitLevel(NewThreshold);

			Node *Left=node[i]->GetLeftNode();
			Node *Right=node[i]->GetRightNode();

			std::vector <std::vector <int> > NewLeftLeaves, NewRightLeaves;
			std::vector <int> NewLeftSubjectList, NewRightSubjectList;

			for (int j=0;j<Leaves.size();j++)
			{
			  if (obs->GetX(Leaves[j][0], NewVariable) <= NewThreshold)
			    NewLeftLeaves.push_back(Leaves[j]);
			  if (obs->GetX(Leaves[j][0], NewVariable) > NewThreshold)
			    NewRightLeaves.push_back(Leaves[j]);
			};

			for (int j=0;j<SubjectList.size();j++)
			{
			  if (obs->GetX(SubjectList[j], NewVariable) <= NewThreshold)
			    NewLeftSubjectList.push_back(SubjectList[j]);
			  if (obs->GetX(SubjectList[j], NewVariable) > NewThreshold)
			    NewRightSubjectList.push_back(SubjectList[j]);
			};

			node.push_back(Left);
			node.push_back(Right);

			LeavesList.push_back(NewLeftLeaves);
			LeavesList.push_back(NewRightLeaves);
		};
	};

	potential= potential + pot;
	Delta *delta = new DeltaBasic(NewTree);
	return delta;
};

void ProposalBasicRadical::ExtendSplittingVariables(std::vector <int> CurrentSplittingVariables, std::vector <int> &NewSplittingVariables, double &pot, bool extend) const
{
	if (extend==true)
	{
		pot+=0;
		NewSplittingVariables=CurrentSplittingVariables;
		return;
	}
	else
	{
		pot+=0;
		return;
	};
};

void ProposalBasicRadical::Sort(std::vector <double> &max, std::vector <double> &min) const
{
	// Sort max and min accoring to the values in max;
	// Currently using Bubble algorithm. Could use more efficient one;
	int i, j;
	if (max.size()!=min.size())
	{
		Rcout<<"Something wrong w ith Sort. Terminated!\n";
	};
	for (i=0;i<max.size()-1;i++)
	{
		double temp=max[i];
		int p=i;
		for (j=i+1;j<max.size();j++)
		{
			if (max[j]<temp) // Sort ascendently
			{
				temp=max[j];
				p=j;
			};
		};
		temp=max[p];
		max[p]=max[i];
		max[i]=temp;

		temp=min[p];
		min[p]=min[i];
		min[i]=temp;
	};
	return;
};

void ProposalBasicRadical::ShowLeaves(std::vector <std::vector <int> > Leaves) const
{
	int i,j;
	for (i=0;i<Leaves.size();i++)
	{
		for (j=0;j<Leaves[i].size();j++)
			Rcout<<Leaves[i][j]<<"\t";
		Rcout<<endl;
	};
	return;
};

void ProposalBasicRadical::ShowSplittingVariables(std::vector <int> SplittingVariables) const
{
	int i;
	for (i=0;i<SplittingVariables.size();i++)
		Rcout<<SplittingVariables[i]<<"\t";
	Rcout<<endl;
};

void ProposalBasicRadical::ShowSubjectList(std::vector <int> SubjectList) const
{
	int i;
	for (i=0;i<SubjectList.size();i++)
		Rcout<<SubjectList[i]<<"\t";
	Rcout<<endl;
};

void ProposalBasicRadical::ShowMaxMin(std::vector <double> Max, std::vector <double> Min) const
{
	int i;
	for (i=0;i<Max.size();i++)
		Rcout<<Max[i]<<"\t"<<Min[i]<<endl;
	Rcout<<endl;
};
