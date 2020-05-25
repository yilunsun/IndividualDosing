#ifndef NODETREE_H
#define NODETREE_H

#ifndef _DIAGRAMS_NODETREE_H
#define _DIAGRAMS_NODETREE_H
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <string>


#include "Node.h"
#include "Random.h"
#include "Density.h"


using namespace std;
using namespace Rcpp;




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
  int GetMiniNodeSize() const;
  std::vector <Node *> GetLeaves();
  List GetVal(); //get value of a tree and optimal label for each node
  std::vector <std::vector <double> > GetLeavesObservations();
  std::vector <std::vector <int> > GetLeavesSubjects();
  std::vector<double> GAM(const NumericVector& a, const NumericVector& ps, const NumericVector& y, const NumericVector& a_out);
  NumericVector GCBPS(const NumericMatrix& x, const NumericVector& y);
  
  void vector_max(std::vector<double> v, double &max, int &imax);
  void checkTree(void) const;
  
 private:
  double potential;
    
  int nLeaf;
  int nLevel;
};

#endif
#endif
