#ifndef DP_H
#define DP_H

#ifndef _DIAGRAMS_DP_H
#define _DIAGRAMS_DP_H

//#include "DeltaGrow.h"

class ModelTreeStructurePinBall;
class ModelSplitVariableIndependent;
class ModelLikelihoodMultinomial;

class DP
{
public:
	static double PotentialDifference(const NodeTree &oldTree,const Delta &delta,
		const ModelTreeStructurePinBall &model,
		Random &ran);

	static double PotentialDifference(const NodeTree &oldTree,const Delta &delta,
		const ModelSplitVariableIndependent &model,
		Random &ran);

  static double PotentialDifference(const NodeTree &oldTree,const Delta &newTree,
    const ModelLikelihoodMultinomial &model,
    Random &ran);

protected:
};


#endif
#endif
