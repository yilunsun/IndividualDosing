#ifndef MODEL_H
#define MODEL_H

#ifndef _DIAGRAMS_MODEL_H
#define _DIAGRAMS_MODEL_H

#include <vector>
#include "Random.h"

using namespace std;

class NodeTree;

class Delta;

class Model
{
 public:
  Model(void);
  virtual ~Model(void);
  virtual double Potential(const NodeTree &tree,Random &ran) const = 0;
  virtual double PotentialDifference(const NodeTree &oldTree,
	  const Delta &newTree,
	  Random &ran) const = 0;
  virtual Model *Copy(void) const = 0;

 private:

};



inline Model::Model(void)
{
  return;
};


inline Model::~Model(void)
{
  return;
};



#endif
#endif
