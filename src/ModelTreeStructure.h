#ifndef MODELTREESTRUCTURE_H
#define MODELTREESTRUCTURE_H

#ifndef _DIAGRAMS_MODELTREESTRUCTURE_H
#define _DIAGRAMS_MODELTREESTRUCTURE_H

#include "Model.h"

class ModelTreeStructure : public Model
{
 public:
  ModelTreeStructure(void);
  virtual ~ModelTreeStructure(void);
  virtual double Potential(const NodeTree &tree,Random &ran) const = 0;
  virtual double PotentialDifference(const NodeTree &oldTree,
				     const Delta &newTree,Random &ran) const = 0;
  virtual NodeTree *Simulate(Random &ran, Observation *obs) const = 0;
  virtual Model *Copy(void) const = 0;

 private:

};



inline ModelTreeStructure::ModelTreeStructure(void) : Model()
{
  return;
};


inline ModelTreeStructure::~ModelTreeStructure(void)
{
  return;
};


#endif
#endif
