#ifndef PROPOSALSWAP_H
#define PROPOSALSWAP_H

#ifndef _DIAGRAMS_PROPOSALSWAP_H
#define _DIAGRAMS_PROPOSALSWAP_H

#include "Proposal.h"


class ProposalSwap : public Proposal 
{
 public:
  ProposalSwap();
  // prob and dens are for selecting new splitting variable and rule after the tree grows.
  ~ProposalSwap(void);

  Delta *ProposeChange(const NodeTree &tree,double &potential,Random &ran) const;
  // should return a DeltaPrune or a DeltaGrow

 private:

};


#endif
#endif