#ifndef MODELLIKELIHOODMULTINOMIAL_H
#define MODELLIKELIHOODMULTINOMIAL_H

#ifndef _DIAGRAMS_MODELMULTINOMIAL_H
#define _DIAGRAMS_MODELMULTINOMIAL_H

#include "ModelLikelihood.h"
#include "Observation.h"

class Node;

class ModelLikelihoodMultinomial : public ModelLikelihood
{
public:
    ModelLikelihoodMultinomial(std::vector<double> alpha, int k);
    ModelLikelihoodMultinomial();
    ~ModelLikelihoodMultinomial(void);

    double PotentialDifference(const NodeTree &oldTree,const Delta &newTree,
                               Random &ran) const;
    Model *Copy(void) const;
    double Potential(const NodeTree &tree,Random &ran) const;
    double PotentialSubTree(const Node &tree,Random &ran) const;

    double PotentialDifferencePrune(const Node *newLeaf,Random &ran) const;
    double PotentialDifferenceGrow(const Node *oldLeaf,int newVar,double newTau,Random &ran) const;

    double PotentialOneLeaf(const std::vector<int> &subject,const Observation *obs,Random &ran) const;
    double IntApprox(std::vector <double> Y, Random &ran) const;
    int GetK(void){return k;}
private:
    std::vector<double> alpha;
    int k;
};


#endif
#endif
