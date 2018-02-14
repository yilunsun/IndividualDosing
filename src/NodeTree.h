#ifndef NODETREE_H
#define NODETREE_H

#ifndef _DIAGRAMS_NODETREE_H
#define _DIAGRAMS_NODETREE_H

#include <iostream>
#include <string>

using namespace std;
using namespace Rcpp;
#include <Rcpp.h>
#include "Node.h"
#include "Random.h"
#include "Density.h"

class Model;

class NodeTree : public Node
{
 public:
  NodeTree();
  NodeTree(const Observation *obs);
  NodeTree(Node &root);
  ~NodeTree(void);

  using Node::GetParentNode;
  Node *GetParentNode(void) {return NULL;};
  void SetParentNode(Node *node);    // a call to this function should
                                     // give an error message
  double Pot(const Model &modelStructure, const Model &modelVar,  const Model &modelLikelihood, Random &ran);
  double LogLikelihood(const Model &modelLikelihood, Random &ran);
  NodeTree *CopyTree(void) const;
  double GetPot() {return potential;};
  bool Misclassfication(int &LeavesSize, int &Error, Density *like) const;
  std::vector <int> GetSplittingVariableSet() const;

  int GetSize() const;
  std::vector <Node *> GetLeaves();
  std::vector <std::vector <double> > GetLeavesObservations();
  std::vector <std::vector <int> > GetLeavesSubjects();

 private:
  double potential;

  int nLeaf;
  int nLevel;
};

#endif
#endif