#ifndef MODELSPLITVARIABLE_H
#define MODELSPLITVARIABLE_H

#ifndef _DIAGRAMS_MODELSPLITVARIABLE_H
#define _DIAGRAMS_MODELSPLITVARIABLE_H

#include "Model.h"

class ModelSplitVariable : public Model
{
 public:
  ModelSplitVariable(void);
  virtual ~ModelSplitVariable(void);
  virtual double Potential(const NodeTree &tree,Random &ran) const = 0;
  virtual double PotentialDifference(const NodeTree &oldTree,
				     const Delta &newTree,Random &ran) const = 0;
  virtual void Simulate(NodeTree &tree,Random &ran) const = 0;
  virtual Model *Copy(void) const = 0;

 private:

};



inline ModelSplitVariable::ModelSplitVariable(void)
{
  return;
};


inline ModelSplitVariable::~ModelSplitVariable(void)
{
  return;
};

#endif
#endif