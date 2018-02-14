#include <string>
#include <algorithm>
#include "Random.h"
#include "Observation.h"
#include "Model.h"

#include "ModelLikelihoodMultinomial.h"
#include "Density.h"
#include "NodeTree.h"
#include "Diagnose.h"


Diagnose::Diagnose()
{
  return;
};

Diagnose::Diagnose(std::vector <NodeTree *> Sample0, Density &density, Model &priorTreeStructure, Model &priorSplitVariable, Model &likelihood)
{
	like = density.Copy();
	mPriorStructure = priorTreeStructure.Copy();
	mPriorSplit = priorSplitVariable.Copy();
	mLike = likelihood.Copy();

	Sample=Sample0;
};


std::vector <NodeTree *> Diagnose::Choose(std::vector <NodeTree *> Sample0, Random &ran)
{
  std::vector <double> weight(Sample0.size(), 0);
  std::vector <NodeTree *> Temp;
  std::vector <int> index;
	for (int i=0;i<weight.size();i++)
	{
		weight[i]=-Sample[i]->Pot(*mPriorStructure, *mPriorSplit, *mLike, ran);
		index.push_back(i);
	};

	for (int i=0;i<weight.size()-1;i++)
	{
		double max=weight[i]; int pointer=i;
		int temp=0;

		for (int j=i+1;j<weight.size();j++)
			if (weight[j]>max)
			{
				max=weight[j]; pointer=j;
			};

		weight[pointer]=weight[i]; weight[i]=max;
		temp=index[i]; index[i]=index[pointer]; index[pointer]=temp;
	};

	int size=(int) Sample0.size()*0.1;
	for (int i=0;i<size;i++)
		Temp.push_back(Sample0[i]);
	return Temp;

};


Diagnose::~Diagnose()
{
	delete like;
	delete mPriorStructure;
	delete mPriorSplit;
	delete mLike;
};


List Diagnose::Scoring(Random &ran)
{
  std::vector <int> Observees;

  // Sort tree first
  this->Sort(ran);

  // Select observees.
  // Modified on Dec 11th. To add only a few observees. Changed from Sample[0]->GetObservation()->GetN() to 4

  for (int i=0;i<Sample[0]->GetObservation()->GetN();i++)
    Observees.push_back(i);

  std::vector <std::vector <double> > llkhd = LogLikelihood(ran);
  std::vector <double> ppot=PosteriorPotential(ran);
  std::vector <int> treesize;

  for(int i=0; i<Sample.size();i++)
    treesize.push_back(Sample[i]->GetSize());

  int which_max=std::distance(ppot.begin(),std::min_element(ppot.begin(),ppot.end()));

  List MAPtree = Sample[which_max]->DumpDotFile();

  return List::create(Named("loglikelihood")=llkhd,
                      Named("posteriorpotential")=ppot,
                      Named("treesize")=treesize,
                      Named("MAPtree")=MAPtree);
};


std::vector <std::vector <double> > Diagnose::SumofSquare()
{
  std::vector <std::vector <double> > LeavesObs;

  std::vector <std::vector <double> > SSquare;

	for (int i=0;i<Sample.size();i++)
	{
		double Sum2=0;
		LeavesObs=Sample[i]->GetLeavesObservations();
		for (int j=0;j<LeavesObs.size();j++)
		{
			double sum=0, mean;
			for (int k=0;k<LeavesObs[j].size();k++)
				sum=sum+LeavesObs[j][k];
			mean=sum/((double) LeavesObs[j].size());

			sum=0;
			for (int k=0;k<LeavesObs[j].size();k++)
				sum=sum+(LeavesObs[j][k]-mean)*(LeavesObs[j][k]-mean);
			Sum2=Sum2+sum;
		};

		std::vector <double> temp;
		for (int k=0;k<2;k++){
		  temp.push_back(Sample[i]->GetSize());
		  temp.push_back(Sum2);
		}
		SSquare.push_back(temp);
	};

	return SSquare;
};


std::vector <double> Diagnose::PosteriorPotential(Random &ran)
{
  std::vector <double> pot;

	for (int i=0;i<Sample.size();i++)
	{
	  double tempresult=Sample[i]->Pot(*mPriorStructure, *mPriorSplit, *mLike, ran);
	  pot.push_back(tempresult);
	}

	return pot;
};


std::vector <std::vector <int> > Diagnose::ImportantCovariates(Random &ran)
{
	int max=0;

	int p=Sample[0]->GetObservation()->GetP();

	std::vector <int> temp(p,0);

	std::vector <std::vector <int> > Weight(Sample[0]->GetObservation()->GetN()+1, temp);

	for (int i=0;i<Sample.size();i++)
	{
	  std::vector <int> Splitting=Sample[i]->GetSplittingVariableSet();
	  std::vector <int> flag(Sample[0]->GetObservation()->GetP(),0);

		int size=Sample[i]->GetSize();
		if (size>max)
			max=size;
		for (int j=0;j<Splitting.size();j++)
			flag[Splitting[j]]=1;

		for (int j=0;j<flag.size();j++)
			if (flag[j]==1)
				Weight[size][j]++;
	};

	return Weight;
};

//Modified on Dec 15
List Diagnose::ImportantCovariates2(Random &ran)
{
	int max=0;

	int p=Sample[0]->GetObservation()->GetP();
	std::vector <int> temp(Sample[0]->GetObservation()->GetP(),0);
	std::vector <std::vector <int> > Weight(Sample[0]->GetObservation()->GetN()+1, temp);
	std::vector <std::vector <int> > Table(p, temp);
	std::vector <std::vector <int> > Count;

	for (int i=0;i<Sample.size();i++)
	{
	  std::vector <int> Splitting=Sample[i]->GetSplittingVariableSet();
	  std::vector <int> flag(Sample[0]->GetObservation()->GetP(),0);

		int size=Sample[i]->GetSize();
		if (size>max)
			max=size;
		for (int j=0;j<Splitting.size();j++)
			flag[Splitting[j]]++;

		for (int j=0;j<p;j++){
		  std::vector <int> temp(3);
		  for (int k=0;k<3;k++){
		    temp[0] = i;
		    temp[1] = j;
		    temp[2] = flag[j];
		  }
		  Count.push_back(temp);
		}

		for (int j=0;j<flag.size();j++)
			if (flag[j]!=0)
				Weight[size][j]++;

		for (int j=0;j<p;j++)
			for (int k=0;k<p;k++)
				if (flag[j]!=0 && flag[k]!=0)
					Table[j][k]++;
	};

	return List::create(Named("Covariates")=Weight,
                      Named("ColTable")=Table,
                      Named("Count")=Count);
};


std::vector <int> Diagnose::GetSizeofTree()
{
  std::vector <int> Result;
	for (int i=0;i<Sample.size();i++)
	{
		Result.push_back(Sample[i]->GetSize());
	};
	return Result;
};


int Diagnose::Seek(std::vector <std::vector <int> > &SubjectList, int j)
{
	int Result;
	bool flag=false;
	for (int i=0; i<SubjectList.size();i++)
		for (int k=0; k<SubjectList[i].size(); k++)
			if (SubjectList[i][k]==j)
			{
				if (! flag)
				{
					Result=i;
					flag=true;
				}
				else
				{
					Rcout<<"Duplicate Subject"<<endl;
				};
			};
	if (! flag)
	{
		Rcout<<"Subject not found.\n";
	};
	return Result;
};

Node * Diagnose::FindLeaf(Node * node, std::vector <double> &x)
{
	while (node->GetLeftNode()!=NULL)
	{
		if(node->GetSplitVariable()>=x.size() || node->GetSplitVariable()<0)
		  Rcout<<"Something wrong about splitting variable in diagnose.cpp!"<<endl;

		if (x[node->GetSplitVariable()] <= node->GetSplitLevel() )
			node=node->GetLeftNode();
		else
			node=node->GetRightNode();
	};
	return node;
};

std::vector <std::vector <double> > Diagnose::LogLikelihood(Random &ran)
{
  std::vector <std::vector <double> > llikelihood;

	for (int i=0;i<Sample.size();i++)
	{
	  std::vector <double> temp;
	  temp.push_back(i+1);
	  temp.push_back(Sample[i]->GetSize());
	  temp.push_back(Sample[i]->LogLikelihood(*mLike, ran));
	  llikelihood.push_back(temp);
	};

	return llikelihood;
};

bool Diagnose::Sort(Random &ran)
{
	//Use potentila or likelihood as the criterion to sort the tree.
	for (int i=0;i<Sample.size();i++)
	{
		//double temp=Sample[i]->Pot(*mPriorStructure, *mPriorSplit, *mLike, ran);
		double temp=-Sample[i]->LogLikelihood(*mLike, ran);
		Potential.push_back(temp);
	};

	for (int i=0;i<Potential.size();i++)
	{
		double min=Potential[i]; // Pot is negative log probability. If log likelihood is used, we use negative log likelihood here.
		int index=i;
		for (int j=i+1;j<Potential.size();j++)
		{
			if (Potential[j]<min)
			{
				min=Potential[j];
				index=j;
			};
		};
		if (min<Potential[i])
		{
			NodeTree *Temp=Sample[i];
			Sample[i]=Sample[index];
			Sample[index]=Temp;
			double temp_pot=Potential[i];
			Potential[i]=Potential[index];
			Potential[index]=temp_pot;
		};
	};
	return true;
};

