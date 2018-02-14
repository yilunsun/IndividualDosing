#ifndef PROPOSALCHANGE_H
#define PROPOSALCHANGE_H

#ifndef _DIAGRAMS_PROPOSALCHANGE_H
#define _DIAGRAMS_PROPOSALCHANGE_H

#include "Proposal.h"
#include "Density.h"
#include "Delta.h"
#include "NodeTree.h"


class ProposalChange : public Proposal
{
 public:
  ProposalChange(const std::vector<double> &probability,
		 const std::vector<Density *> &density);
  ~ProposalChange(void);

  Delta *ProposeChange(const NodeTree &tree,double &potential,Random &ran) const;
  // should return a DeltaChange

 private:
  std::vector<double> prob;
  std::vector<Density *> dens;
};


#endif
#endif
