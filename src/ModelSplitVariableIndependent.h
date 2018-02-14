#ifndef MODELSPLITVARIABLEINDEPENDENT_H
#define MODELSPLITVARIABLEINDEPENDENT_H

#ifndef _DIAGRAMS_MODELSPLITVARIABLEINDEPENDENT_H
#define _DIAGRAMS_MODELSPLITVARIABLEINDEPENDENT_H

#include "ModelSplitVariable.h"
#include "Density.h"


class ModelSplitVariableIndependent : public ModelSplitVariable
{
 public:
  ModelSplitVariableIndependent(const std::vector<double> &probability,
				const std::vector<Density *> &density);
  ~ModelSplitVariableIndependent(void);

  double PotentialDifference(const NodeTree &oldTree,
			     const Delta &newTree,Random &ran) const;
  void Simulate(NodeTree &tree,Random &ran) const;
  Model *Copy(void) const;
  double Potential(const NodeTree &tree,Random &ran) const;
  double PotentialOneSplit(int nr,double tau,Random &ran) const;

 private:
  std::vector<double> prob;
  std::vector<Density *> dens;
};

#endif
#endif
