#ifndef MODELLIKELIHOOD_H
#define MODELLIKELIHOOD_H

#ifndef _DIAGRAMS_MODELLIKELIHOOD_H
#define _DIAGRAMS_MODELLIKELIHOOD_H

#include "Model.h"

class ModelLikelihood : public Model
{
 public:
  ModelLikelihood(void);
  ~ModelLikelihood(void);
  double Potential(const NodeTree &tree,Random &ran) const {return 0;};
  double PotentialDifference(const NodeTree &oldTree,
				     const Delta &newTree,Random &ran) const {return 0;};
  Model *Copy(void) const {return NULL;};



 private:

};



inline ModelLikelihood::ModelLikelihood(void)
{
  return;
};



inline ModelLikelihood::~ModelLikelihood(void)
{
  return;
};

#endif
#endif