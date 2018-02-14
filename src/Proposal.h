#ifndef PROPOSAL_H
#define PROPOSAL_H

#ifndef _DIAGRAMS_PROPOSAL_H
#define _DIAGRAMS_PROPOSAL_H



class Proposal
{
 public:
  Proposal(void);
  virtual ~Proposal(void);

  virtual Delta *ProposeChange(const NodeTree &tree,double &potential,
			       Random &ran) const = 0;
  
 private:

};


inline Proposal::Proposal(void)
{
  return;
};



inline Proposal::~Proposal(void)
{
  return;
};


#endif
#endif