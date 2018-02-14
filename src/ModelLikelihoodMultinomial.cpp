#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "NodeTree.h"
#include "Delta.h"
#include "DP.h"
#include "ModelLikelihoodMultinomial.h"

#include <numeric>

ModelLikelihoodMultinomial::ModelLikelihoodMultinomial(std::vector<double> alpha, int k)
{
    //
    // the prior is:  p(k-1 dimensional) ~ Dirichlet(alpha,length(alpha)=k)
    // and alpha is (alpha_0, alpha_1,...), in beta, it is (beta, alpha)

    this->alpha = alpha;
    this->k = k;

    return;
};

ModelLikelihoodMultinomial::ModelLikelihoodMultinomial()
{
    return;
};



ModelLikelihoodMultinomial::~ModelLikelihoodMultinomial(void)
{
    return;
};



double ModelLikelihoodMultinomial::PotentialDifference(const NodeTree &oldTree,
                                                    const Delta &newTree,
                                                    Random &ran) const
{
    double pot = DP::PotentialDifference(oldTree,newTree,*this,ran);

    return pot;
};




Model *ModelLikelihoodMultinomial::Copy(void) const
{
    Model *model = new ModelLikelihoodMultinomial(alpha, k);

    return model;
};

double ModelLikelihoodMultinomial::IntApprox(std::vector <double> Y, Random &ran) const
{
    int l=Y.size();

    std::vector <int> s1;

    int s;

    double pot=0;

    int i;

    for (i=0; i<k; i++)
    {
        int count=0;

        for (s=0; s<l; s++)
        {
            if (Y[s]==i) count++;
        };

        s1.push_back(count);
    }

    if(std::accumulate(s1.begin(), s1.end(), 0.0)!=l)
      Rcout<<"IntApprox: something wrong with the likelihood!"<<endl;

    for (i=0;i<k;i++)
    {
        for (s=0;s<s1[i];s++)
        {
            pot += -log(alpha[i] + s);
        }
    }

    for (s=0;s<l;s++)
        pot -= -log(std::accumulate(alpha.begin(), alpha.end(), 0.0)+s);
    return pot;
};


double ModelLikelihoodMultinomial::Potential(const NodeTree &tree,Random &ran) const
{
    double pot = 0.0;

    std::vector<const Node *> node;

    node.push_back(&tree);
    int i;
    for (i = 0; i < node.size(); i++)
    {
        if (node[i]->GetRightNode() != NULL)  // interior node
        {
            node.push_back(node[i]->GetRightNode());
            node.push_back(node[i]->GetLeftNode());
        }
        else  // a leaf, contributes to the potential
        {
            std::vector<int> subject(node[i]->GetSubjectList());
            const Observation *obs = node[i]->GetObservation();

            std::vector <int> s1;

            int l = subject.size();

            int s;

            std::vector <double> Y;

            int j;

            for (j=0; j<k; j++)
            {
                int count = 0;

                for (s = 0; s < l; s++)
                {
                    int nr = subject[s];
                    Y.push_back(obs->GetY(nr));
                    if (obs->GetY(nr)==j) count++;
                }

                s1.push_back(count);
            }

            if(std::accumulate(s1.begin(), s1.end(), 0.0)!=l)
              Rcout<<"Potential: something wrong with the likelihood!"<<endl;

            //double pottemp=IntApprox(Y, ran);

            for (j=0;j<k;j++)
            {
                for (s=0;s<s1[j];s++)
                {
                    pot += -log(alpha[j] + s);
                }
            }

            for (s=0;s<l;s++)
                pot -= -log(std::accumulate(alpha.begin(), alpha.end(), 0.0)+s);

        }
    }

    return pot;
};







double ModelLikelihoodMultinomial::PotentialSubTree(const Node &subtree,Random &ran) const
{
    double pot = 0.0;

    std::vector<const Node *> node;

    node.push_back(&subtree);
    int i;
    for (i = 0; i < node.size(); i++)
    {
        if (node[i]->GetRightNode() != NULL)  // interior node
        {
            node.push_back(node[i]->GetRightNode());
            node.push_back(node[i]->GetLeftNode());
        }
        else  // a leaf, contributes to the potential
        {
            std::vector<int> subject(node[i]->GetSubjectList());
            const Observation *obs = node[i]->GetObservation();

            pot += PotentialOneLeaf(subject,obs,ran);
        }
    }

    return pot;
};



double ModelLikelihoodMultinomial::PotentialOneLeaf(const std::vector<int> &subject,const Observation *obs,Random &ran) const
{
    std::vector <int> s1;

    int l = subject.size();

    int s;

    int i;

    for (i=0; i<k; i++)
    {
        int count=0;

        for (s = 0; s < l; s++)
        {
            int nr = subject[s];
            if (obs->GetY(nr)==i) count++;
        }

        s1.push_back(count);
    }

    if (std::accumulate(s1.begin(), s1.end(), 0.0)!=l){
      Rcout<<"PotentialOneLeaf: something wrong with the likelihood!"<<endl;
    }

    double pot = 0.0;

    for (i=0;i<k;i++)
    {
        for (s=0;s<s1[i];s++)
        {
            pot += -log(alpha[i] + s);
        }
    }

    for (s=0;s<l;s++)
        pot -= -log(std::accumulate(alpha.begin(), alpha.end(), 0.0)+s);

    return pot;
};



double ModelLikelihoodMultinomial::PotentialDifferencePrune(const Node *newLeaf,Random &ran) const
{
    double pot = 0.0;

    std::vector<int> subjectL(newLeaf->GetLeftNode()->GetSubjectList());
    std::vector<int> subjectR(newLeaf->GetRightNode()->GetSubjectList());
    const Observation *obs = newLeaf->GetLeftNode()->GetObservation();

    pot -= PotentialOneLeaf(subjectL,obs,ran);
    pot -= PotentialOneLeaf(subjectR,obs,ran);

    std::vector<int> subject;
    int i;
    for (i = 0; i < subjectL.size(); i++)
        subject.push_back(subjectL[i]);
    for (i = 0; i < subjectR.size(); i++)
        subject.push_back(subjectR[i]);

    pot += PotentialOneLeaf(subject,obs,ran);

    return pot;
};




double ModelLikelihoodMultinomial::PotentialDifferenceGrow(const Node *oldLeaf,int newVar,double newTau,Random &ran) const
{
    double pot = 0.0;

    std::vector<int> subject(oldLeaf->GetSubjectList());
    const Observation *obs = oldLeaf->GetObservation();

    pot -= PotentialOneLeaf(subject,obs,ran);



    Node *newInterior = oldLeaf->Copy();
    Node *newLeft = new Node;
    newLeft->SetObservation(obs);
    newLeft->SetSplitVariable(-1);
    newLeft->SetMinimumLeafSize(oldLeaf->GetMinimumLeafSize());

    Node *newRight = new Node;
    newRight->SetObservation(obs);
    newRight->SetSplitVariable(-1);
    newRight->SetMinimumLeafSize(oldLeaf->GetMinimumLeafSize());

    newInterior->SetLeftNode(newLeft);  newInterior->SetRightNode(newRight);
    newInterior->SetSplitVariable(newVar);
    newInterior->SetSplitLevel(newTau);
    if (!newInterior->UpdateSubjectList())
        Rcout << "debug\n";    // should never happen!

    pot += PotentialOneLeaf(newLeft->GetSubjectList(),obs,ran);
    pot += PotentialOneLeaf(newRight->GetSubjectList(),obs,ran);

    newInterior->RemoveSubTree();
    delete newInterior;
    newInterior = NULL;

    return pot;
};
