#ifndef DIAGNOSE_H
#define DIAGNOSE_H

#ifndef _DIAGRAMS_DIAGNOSE_H
#define _DIAGRAMS_DIAGNOSE_H

#include <iostream>
#include <vector>
#include <Rcpp.h>
#include "Random.h"
#include "NodeTree.h"
#include "Density.h"
using namespace std;
using namespace Rcpp;

class Diagnose
{
public:
	Diagnose(void);
	Diagnose(std::vector <NodeTree *> Sample0, Density &density, Model &priorTreeStructure,Model &priorSplitVariable,Model &likelihood);
	std::vector <int> GetSizeofTree();
	virtual ~Diagnose(void);

	List Scoring(Random &ran);
	std::vector <double> PosteriorPotential(Random &ran);
	std::vector <std::vector <double> > LogLikelihood(Random &ran);
	std::vector <std::vector <int> > ImportantCovariates(Random &ran);
	List ImportantCovariates2(Random &ran);
	std::vector <std::vector <double> > CV(Observation *Test, Random &ran);
	std::vector <std::vector <double> > ErrorTable();
	std::vector <std::vector <double> > SumofSquare();
	std::vector <NodeTree *> Choose(std::vector <NodeTree *> Sample0, Random &ran);

private:
  std::vector <NodeTree *> Sample;
	int Seek(std::vector <std::vector <int> > &SubjectList, int j);
	bool Sort(Random &ran);

	std::vector <double> Potential;

	Node * FindLeaf(Node * node, std::vector <double> &x);
	Density *like;
	Model *mPriorStructure;
	Model *mPriorSplit;
	Model *mLike;
};

#endif
#endif
