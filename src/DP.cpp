#include "Random.h"

#include "Density.h"
#include "Observation.h"
#include "ModelTreeStructurePinBall.h"
#include "ModelSplitVariableIndependent.h"
#include "ModelLikelihoodMultinomial.h"
#include "NodeTree.h"

#include "DeltaBasic.h"
#include "DeltaChange.h"
#include "DeltaPrune.h"
#include "DeltaGrow.h"
#include "DP.h"

double DP::PotentialDifference(const NodeTree &oldTree,const Delta &newTree,
							   const ModelTreeStructurePinBall &model,
							   Random &ran)
{
	int type = newTree.GetType();

	if (type == DELTA_BASIC)
	{
		const DeltaBasic *delta = (const DeltaBasic *) &newTree;
		const NodeTree *nn = delta->GetNewTree();

		double pot = model.Potential(*nn,ran);
		pot -= model.Potential(oldTree,ran);

		return pot;
	}
	else if (type == DELTA_CHANGE)
	{
		return 0.0;
	}
	if (type == DELTA_PRUNE)
	{
		const DeltaPrune *delta = (const DeltaPrune *) &newTree;
		const Node *newLeaf = delta->GetNewLeaf();
		double pot = model.PotentialDifferencePrune(oldTree,newLeaf,ran);

		return pot;
	}
	if (type == DELTA_GROW)
	{
		const DeltaGrow *delta = (const DeltaGrow *) &newTree;
		const Node *oldLeaf = delta->GetOldLeaf();
		double pot = model.PotentialDifferenceGrow(oldTree,oldLeaf,ran);

		return pot;
	};

	Rcout << "ERROR: Unknown Delta-type " << type << " found in function \"DP::PotentialDifference\"\n";

	return 0.0;
};





double DP::PotentialDifference(const NodeTree &oldTree,const Delta &newTree,
							   const ModelSplitVariableIndependent &model,
							   Random &ran)
{
	int type = newTree.GetType();

	if (type == DELTA_BASIC)
	{
		const DeltaBasic *delta = (const DeltaBasic *) &newTree;
		const NodeTree *nn = delta->GetNewTree();

		double pot = model.Potential(*nn,ran);
		pot -= model.Potential(oldTree,ran);

		return pot;
	}
	else if (type == DELTA_CHANGE)
	{
		const DeltaChange *delta = (const DeltaChange *) &newTree;

		double pot = model.PotentialOneSplit(delta->GetVarNew(),delta->GetTauNew(),ran);
		pot -= model.PotentialOneSplit(delta->GetVarOld(),delta->GetTauOld(),ran);

		return pot;
	}
	if (type == DELTA_PRUNE)
	{
		const DeltaPrune *delta = (const DeltaPrune *) &newTree;
		const Node *newLeaf = delta->GetNewLeaf();
		int oldVar = newLeaf->GetSplitVariable();
		double oldTau = newLeaf->GetSplitLevel();

		double pot = - model.PotentialOneSplit(oldVar,oldTau,ran);

		return pot;
	}
	if (type == DELTA_GROW)
	{
		const DeltaGrow *delta = (const DeltaGrow *) &newTree;
		int newVar = delta->GetNewVar();
		double newTau = delta->GetNewTau();

		double pot = model.PotentialOneSplit(newVar,newTau,ran);

		return pot;
	};

	Rcout << "ERROR: Unknown Delta-type " << type << " found in function \"DP::PotentialDifference\"\n";

	return 0.0;
};

double DP::PotentialDifference(const NodeTree &oldTree,const Delta &newTree,
                               const ModelLikelihoodMultinomial &model,
                               Random &ran)
{
    int type = newTree.GetType();

    if (type == DELTA_BASIC)
    {
        const DeltaBasic *delta = (const DeltaBasic *) &newTree;
        const NodeTree *nn = delta->GetNewTree();

        double pot = model.Potential(*nn,ran);
        pot -= model.Potential(oldTree,ran);

        return pot;
    }
    else if (type == DELTA_CHANGE)
    {
        const DeltaChange *delta = (const DeltaChange *) &newTree;

        double pot = model.PotentialSubTree(*(delta->GetNodeNew()),ran);
        pot -= model.PotentialSubTree(*(delta->GetNodeOld()),ran);

        return pot;
    }
    if (type == DELTA_PRUNE)
    {
        const DeltaPrune *delta = (const DeltaPrune *) &newTree;
        const Node *newLeaf = delta->GetNewLeaf();

        double pot = model.PotentialDifferencePrune(newLeaf,ran);

        return pot;
    }
    if (type == DELTA_GROW)
    {
        const DeltaGrow *delta = (const DeltaGrow *) &newTree;
        const Node *oldLeaf = delta->GetOldLeaf();
        int newVar = delta->GetNewVar();
        double newTau = delta->GetNewTau();

        double pot = model.PotentialDifferenceGrow(oldLeaf,newVar,newTau,ran);

        return pot;
    };


    Rcout << "ERROR: Unknown Delta-type " << type << " found in function \"DP::PotentialDifference\"\n";

    return 0.0;
};


