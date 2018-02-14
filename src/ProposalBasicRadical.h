#ifndef PROPOSALRADICALBASIC_H
#define PROPOSALRADICALBASIC_H

#ifndef _DIAGRAMS_PROPOSALBASICRADICAL_H
#define _DIAGRAMS_PROPOSALBASICRADICAL_H

#include "Proposal.h"
#include "ChangeLeaves.h"
#include "Delta.h"
#include "NodeTree.h"


class ProposalBasicRadical : public Proposal
{
 public:
    ProposalBasicRadical(const ChangeLeaves &changeLeaves);
    ~ProposalBasicRadical(void);

    Delta *ProposeChange(const NodeTree &tree,double &potential,Random &ran) const;
 protected:
    std::vector <std::vector <int> > GetLeavesConfiguration(const NodeTree &tree, std::vector <int> &CurrentSplittingVariables) const;
    void ExtendSplittingVariables(std::vector <int> CurrentSplittingVariables, std::vector <int> &NewSplittingVariables, double &pot, bool extend) const;
    void Sort(std::vector <double> &max, std::vector <double> &min) const;

    //all test functions
    void ShowLeaves(std::vector <std::vector <int> > Leaves) const;
    void ShowSplittingVariables(std::vector <int> SplittingVariables) const;
    void ShowSubjectList(std::vector <int> SubjectList) const;
    void ShowMaxMin(std::vector <double> Max, std::vector <double> Min) const;

 private:
    ChangeLeaves *changeLeaves;
};


#endif
#endif
