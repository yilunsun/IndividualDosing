#ifndef NODE_H
#define NODE_H

#ifndef _DIAGRAMS_NODE_H
#define _DIAGRAMS_NODE_H

#include <vector>
#include <Rcpp.h>
#include "Observation.h"

using namespace std;
using namespace Rcpp;

class NodeTree;

class Node
{
public:
	Node(void);
	virtual ~Node(void);

	void SetObservation(const Observation* obs);  // set local variable obs for this
	// node and for all descendent nodes
	// I change the reference to pointer.
	Node *GetLeftNode(void) const {return left;};
	Node *GetRightNode(void) const {return right;};
	Node *GetRoot(void);
	const Observation* GetObservation(void) const;

	virtual Node *GetParentNode(void) const {return parent;};
	int GetSplitVariable(void) const {return splitVariable;};
	double GetSplitLevel(void) const {return splitLevel;};
	int GetFirstObservation(void) const;  // Return 0 if no observations
	int GetNextObservation(void) const;   // Return 0 if no more observations
	int GetNumberOfLeaf(void) const;      // Return number of leaves in the subtree with
	// this node as root

	Node *Copy(void) const;               // copy this node only
	Node *CopySubTree(void) const;    // copy this node and its descendents and
	// return this as a tree
	void RemoveSubTree(void);       // Delete the descendents of this node. Not the node itself.

	void SetLeftNode(Node *node);
	void SetRightNode(Node *node);
	virtual void SetParentNode(Node *node);
	void SetSplitVariable(int nr);
	void SetSplitLevel(double value);
	void SetSplittingVar(std::vector <int> var);

	std::vector <int> GetSubjectList(void) const;
	std::vector <int> GetSplittingVar(void) const {return Splitting_var;};
	void ExcSplittingVar(int x);
	void SetSubjectList(std::vector <int> SubjectList);
	void SetMinimumLeafSize(int LeafSize) { MinimumLeafSize=LeafSize; };
	int GetMinimumLeafSize() const { return MinimumLeafSize; };
	bool UpdateSubjectList(void);  // update the subject lists for all
	// descendents of this node

	List DumpDotFile(void);

	std::vector <double> GetSubjectObs();

	void CheckConsistency(void) const; //debugging

protected:
	const Observation *obs;
  std::vector <int> SubjectList;
  std::vector <int> Splitting_var; // variables available for splitting
	int MinimumLeafSize;
	Node *left;
	Node *right;
	Node *parent;   // Equal NULL (only) for root

	int splitVariable;
	double splitLevel; //need to think about this part for current coding
};


#endif
#endif
