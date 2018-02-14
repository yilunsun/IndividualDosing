#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
#include <string>
#include <vector>

#include "Observation.h"
#include "Random.h"
#include "Density.h"
#include "ModelTreeStructurePinBall.h"
#include "ModelSplitVariableIndependent.h"
#include "ModelLikelihoodMultinomial.h"
#include "Delta.h"

#include "ProposalChange.h"
#include "ProposalPruneGrow.h"
#include "ProposalSwap.h"
#include "ProposalBasicRadical.h"


#include "MCMC.h"
#include "NodeTree.h"
#include "DensityDirichlet.h"
#include "DensityUnif.h"
#include "DensityDummy.h"
#include "Diagnose.h"
#include "ChangeLeavesIdentity.h"

using namespace Rcpp;
using namespace std;

//only works on classification problems
//cat_num is number of dummy indicators. They are required to be placed on the left of data matrix, and coded as integer vector 0,1,2,...etc.
//x needs to be standardized to be in [0,1], taking equal sign when there is dummy
//y needs to be coded as an integer vector 0,1,2,...etc.

// [[Rcpp::export]]
List BayesianCART(NumericMatrix x, IntegerVector y, int cat_num, bool standardization, int burnin, int Length, int every, int nChain,
                  double size, double shape, //std::vector <double> importance_weight,
                  String prior_leaf = "uniform", int MinimumLeafSize=1, unsigned int seed=123)
{
  Random ran(seed);

  if (standardization) {
    for (int i=cat_num;i<x.ncol();i++){
      x( _, i) = (x( _, i)-Rcpp::min(x( _, i)))/(Rcpp::max(x( _, i)) - Rcpp::min(x( _, i)));
    }
  }

  List Result;

  Observation obs(x,y,cat_num);
  std::vector <NodeTree *> Sample;

  int p = obs.GetP();
  int cat = obs.GetK();

  //this is used in splitting proposal, default set to equal
  //delete in future
  std::vector <double> importance_weight;
  for (int i=0;i<p;i++)
    importance_weight.push_back(1.0);
  //end of deletion
  if (importance_weight.size()!=p)
    Rcout<<"weight spec is wrong!"<<endl;
  double wt_sum=std::accumulate(importance_weight.begin(), importance_weight.end(), 0.0);
  std::vector <double> prob;
  for (int i=0;i<p;i++)
    prob.push_back(importance_weight[i]/wt_sum);

  //used for splitting rule
  std::vector<Density *> density;

  //categorical variables first
  for (int k = 0; k < cat_num; k++){
    density.push_back(new DensityDummy());
  };

  //now the continuous varibles
  for (int k = cat_num; k < p; k++)
    density.push_back(new DensityUnif(0,1));

  //leaf prior, for now we assume each category has the same prior probability
  std::vector<double> alpha;

  double alpha_each;

  if (prior_leaf == "uniform") {
    alpha_each = 1.0;
  } else if (prior_leaf == "noninform"){
    alpha_each = 0.5;
  };

  for (int i=0;i<cat;i++){
    alpha.push_back(alpha_each);
  };

  ModelTreeStructurePinBall mTreeStructure(size,shape,MinimumLeafSize);
  ModelSplitVariableIndependent mSplitVariable(prob,density);

  ModelLikelihoodMultinomial mLikelihood(alpha,cat);
  DensityDirichlet dDensity(alpha,cat);

  for (int l = 0; l < nChain; l++) // For final data analysis, we just run one long chain.
  {
    Rcout<<"Chain #"<<l+1<<endl;
    NodeTree *tree = mTreeStructure.Simulate(ran, &obs, 1);
    tree->SetMinimumLeafSize(MinimumLeafSize);
    mSplitVariable.Simulate(*tree,ran); // When reading tree from external file, comment this line.
    tree->SetObservation(&obs);

    while (!tree->UpdateSubjectList()){
      mSplitVariable.Simulate(*tree,ran);
    }

    Rcout<<"Initialization Completed."<<endl;

    //
    // Initialise the proposal and MCMC classes
    //

    std::vector<Proposal *> proposal1;
    proposal1.push_back(new ProposalChange(prob,density));
    proposal1.push_back(new ProposalPruneGrow(prob, density));
    proposal1.push_back(new ProposalSwap());

    std::vector<Proposal *> proposal2;
    ChangeLeavesIdentity noChange;
    proposal2.push_back(new ProposalBasicRadical(noChange));

    MCMC mcmc(mTreeStructure,mSplitVariable,mLikelihood);
    std::vector <int> Accept(4, 0);
    //
    // Run the Metropolis-Hastings algorithm
    //

    for (int i = 0; i < Length; i++)
    {

      if (i%50==0)
        Rcout<<"Iteration: "<<i<<endl;

      std::vector<int> nAcc;

      for (int kk = 1; kk <= every; kk++)
      {
        tree = mcmc.Iterate(tree,proposal1,25,nAcc,ran);
        tree = mcmc.Iterate(tree,proposal2,1,nAcc,ran);
        tree = mcmc.Iterate(tree,proposal1,25,nAcc,ran);
      }
      if (i>=burnin)
      {
        NodeTree *temp=tree->CopyTree();
        Sample.push_back(temp);
      };
    };
    Diagnose Diag(Sample, dDensity, mTreeStructure, mSplitVariable, mLikelihood);
    std::string I=std::to_string(l);
    std::string chainname="chain"+I;
    Result[chainname]=Diag.Scoring(ran);

    delete tree;
    for (int k=0;k<proposal1.size();k++)
      delete proposal1[k];
    proposal1.clear();
    for (int k=0;k<proposal2.size();k++)
      delete proposal2[k];
    proposal2.clear();
    for (int k=0;k<Sample.size();k++)
      delete Sample[k];
    Sample.clear();
  };

  return Result;
};

