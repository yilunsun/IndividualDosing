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
  
  arma::colvec repc(const NumericVector& input, const int& ct);
  NumericVector seq_length(const double& start, const double& end, const int& length_out);
  NumericVector dirtyGBM(const NumericMatrix& x, const NumericVector& y, const NumericMatrix& newdf);
  List dirtyDensity(const arma::colvec& x);
  arma::colvec dirtySpline(const arma::colvec& x, const arma::colvec& y, const arma::colvec& newdf);
  arma::colvec w_fn(const double& bw, const arma::colvec a_vals, const arma::colvec& a, const int& n);
  arma::colvec dirtyLocpoly(const double& h, const arma::colvec& out, const arma::colvec& a, const arma::colvec& a_out);
  double risk_fn(const double& h, const arma::colvec& a_vals, const arma::colvec& out, const arma::colvec& a, const int& n);
  std::vector<double> dose_kernel_rf_cpp(const NumericVector& y, const NumericVector& a, const NumericMatrix& x, const NumericVector& a_out, 
                                         double a_min, double a_max, const int n_pts);
  void vector_max(std::vector<double> v, double &max, int &imax);
  void checkTree(void) const;
  
 private:
  double potential;
    
  int nLeaf;
  int nLevel;
};

#endif
#endif
