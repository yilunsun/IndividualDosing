#ifndef MCMC_H
#define MCMC_H


#ifndef _DIAGRAMS_MCMC_H
#define _DIAGRAMS_MCMC_H


class MCMC
{
 public:
  MCMC(Model &priorTreeStructure,Model &priorSplitVariable,Model &likelihood);
  ~MCMC(void);

  NodeTree *Iterate(NodeTree *tree,std::vector<Proposal *> proposal,
		    int numberOfIteration,std::vector<int> &nAccept,
		    Random &ran, double T0, int Length, int i) const;

 private:
  Model *priorStructure;
  Model *priorSplit;
  Model *like;
};

#endif
#endif
