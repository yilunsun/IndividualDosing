#ifndef PROPOSALPRUNEGROW_H
#define PROPOSALPRUNEGROW_H

#ifndef _DIAGRAMS_PROPOSALPRUNEGROW_H
#define _DIAGRAMS_PROPOSALPRUNEGROW_H

#include "Proposal.h"


class ProposalPruneGrow : public Proposal
{
 public:
  ProposalPruneGrow(const std::vector<double> &probability, const std::vector<Density *> &density);
  // prob and dens are for selecting new splitting variable and rule after the tree grows.
  ~ProposalPruneGrow(void);

  Delta *ProposeChange(const NodeTree &tree,double &potential,Random &ran) const;
  // should return a DeltaPrune or a DeltaGrow

 private:
  std::vector<double> prob;
  std::vector<Density *> dens;
};


#endif
#endif
