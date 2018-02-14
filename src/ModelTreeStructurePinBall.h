#ifndef MODELTREESTRUCTUREPINBALL_H
#define MODELTREESTRUCTUREPINBALL_H

#ifndef _DIAGRAMS_MODELTREESTRUCTUREPINBALL_H
#define _DIAGRAMS_MODELTREESTRUCTUREPINBALL_H

#include "ModelTreeStructure.h"

class Node;

class ModelTreeStructurePinBall : public ModelTreeStructure
{
public:
	ModelTreeStructurePinBall(double paramNumberOfLeaves,double paramSplit);
	ModelTreeStructurePinBall(double paramNumberOfLeaves,double paramSplit, int MinimumTreeSize);

	~ModelTreeStructurePinBall(void);

	double PotentialDifference(const NodeTree &oldTree,
		const Delta &newTree,Random &ran) const;
	NodeTree *Simulate(Random &ran, Observation *obs) const;
	NodeTree *Simulate(Random &ran, Observation *obs, int LeafSize) const;
	Model *Copy(void) const;

	double PotentialDifferencePrune(const NodeTree &oldTree,const Node *newLeaf,Random &ran) const;
	double PotentialDifferenceGrow(const NodeTree &oldTree,const Node *oldLeaf,Random &ran) const;

	double Potential(const NodeTree &tree,Random &ran) const;
	double GetParamNumberOfLeaves(void) {return lambda;};
	double GetParamSplit(void) {return prob;};

private:

	double lambda;
	double prob;
	int MinimumSize;
	std::vector <double> GetnLeafProbs(int length, Random &ran) const;
};


#endif
#endif
