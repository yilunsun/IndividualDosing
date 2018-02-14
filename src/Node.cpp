#include <cstdlib>
#include "Observation.h"
#include "Node.h"

// [[Rcpp::plugins("cpp11")]]
Node::Node()
{
  left = NULL;
  right = NULL;
  parent = NULL;
  splitVariable = -1;
  return;
};

Node::~Node()
{
  if (left!=NULL)
  {
    left->RemoveSubTree();
    delete left;
    left=NULL;
  };
  if (right!=NULL)
  {
    right->RemoveSubTree();
    delete right;
    right=NULL;
  };
  return;
};

void Node::SetLeftNode(Node *node)
{
  if (node!=NULL) node->SetParentNode(this);
  left=node;
  return;
};

void Node::SetRightNode(Node *node)
{
  if (node!=NULL) node->SetParentNode(this);
  right=node;
  return;
};

void Node::SetSplitVariable(int nr)
{
  splitVariable=nr;
  return;
};

void Node::SetSplitLevel(double value)
{
  splitLevel=value;
  return;
};

void Node::SetObservation(const Observation* obs)
{
  this->obs=obs;
  return;
};

void Node::SetSplittingVar(std::vector <int> var)
{
  Splitting_var = var;
  return;
};

void Node::ExcSplittingVar(int x)
{
  Splitting_var.erase(std::remove(Splitting_var.begin(), Splitting_var.end(), x), Splitting_var.end());
  return;
};


Node* Node::Copy(void) const
{
  Node* NewNode=new Node;
  NewNode->SetObservation(this->GetObservation());
  NewNode->SetSplitLevel(this->GetSplitLevel());
  NewNode->SetSplitVariable(this->GetSplitVariable());
  NewNode->SetSplittingVar(this->GetSplittingVar());
  NewNode->SetSubjectList(this->GetSubjectList());
  NewNode->SetMinimumLeafSize(this->GetMinimumLeafSize());
  //NewNode->SetParentNode(this->GetParentNode());
  return NewNode;
};

Node* Node::CopySubTree(void) const
{
  Node* NewNode=new Node;
  NewNode->SetObservation(this->GetObservation());
  NewNode->SetSplitLevel(this->GetSplitLevel());
  NewNode->SetSplitVariable(this->GetSplitVariable());
  NewNode->SetSubjectList(this->GetSubjectList());
  NewNode->SetSplittingVar(this->GetSplittingVar());
  NewNode->SetMinimumLeafSize(this->GetMinimumLeafSize());
  if (this->GetLeftNode()!=NULL || this->GetRightNode()!=NULL)
  {
    Node* pLeft=this->GetLeftNode()->CopySubTree();
    Node* pRight=this->GetRightNode()->CopySubTree();
    NewNode->SetLeftNode(pLeft);
    NewNode->SetRightNode(pRight);
  }

  return NewNode;
}

void Node::SetParentNode(Node *node)
{
  parent=node;
  Splitting_var=node->Splitting_var;
  return;
};

const Observation* Node::GetObservation(void) const
{
  return obs;
};

Node* Node::GetRoot(void)
{
  if (!this->parent)
    return this;
  else
    return this->GetParentNode()->GetRoot();
};

std::vector <int> Node::GetSubjectList(void) const
{
  return SubjectList;
};

void Node::SetSubjectList(std::vector <int> SubjectList)
{
  this->SubjectList=SubjectList;
  return;
};

int Node::GetNumberOfLeaf(void) const
{
  //
  // This is a quite inefficent implementation og this function.
  // Perhaps it is better to store some extra information
  // to make this more efficient.
  //

  if (right == NULL)
    return 1;

  //
  // This node is not a leaf!
  //

  std::vector<const Node *> node;
  node.push_back(this);
  int nNode = 0;

  int k;
  for (k = 0; k < node.size(); k++)
  {
    if (node[k]->GetRightNode()->GetRightNode() == NULL)
      nNode++;
    else
      node.push_back(node[k]->GetRightNode());

    if (node[k]->GetLeftNode()->GetRightNode() == NULL)
      nNode++;
    else
      node.push_back(node[k]->GetLeftNode());
  }

  return nNode;
};

bool Node::UpdateSubjectList(void)
{
  std::vector<Node *> node;
  node.push_back(this);

  //if (node.size()==1) return 0; // A single root;
  for (int i = 0; i < node.size(); i++)
  {
    if (node[i]->GetRightNode()!= NULL && node[i]->GetLeftNode()!= NULL)
    {
      int x=node[i]->GetSplitVariable();
      if (std::find(Splitting_var.begin(), Splitting_var.end(), x) == Splitting_var.end())
      {
        Rcout<<i<<"-th node:"<<"\n";
        Rcout<<"Splitting_var:";
        for(int bug=0; bug<Splitting_var.size();bug++)
          Rcout << Splitting_var[bug] << "\t";
        Rcout<<"\n x is:"<<x<<"\n";
      }; //Something wrong.
      double tau=node[i]->GetSplitLevel();
      std::vector <int> Current=node[i]->GetSubjectList();
      std::vector <int> NewLeftSubjectList;
      std::vector <int> NewRightSubjectList;
      for (int j=0;j<Current.size();j++){
        if (node[i]->GetObservation()->GetX(Current[j], x)<=tau)
          NewLeftSubjectList.push_back(Current[j]);
        else NewRightSubjectList.push_back(Current[j]);
      };
      if (NewLeftSubjectList.size()<MinimumLeafSize || NewRightSubjectList.size()<MinimumLeafSize)
        return false; // Not a valid division
      node[i]->GetLeftNode()->SetSubjectList(NewLeftSubjectList);
      node[i]->GetRightNode()->SetSubjectList(NewRightSubjectList);

      node.push_back(node[i]->GetRightNode());
      node.push_back(node[i]->GetLeftNode());
    };
  };

  return true;
};


void Node::RemoveSubTree(void)
{
  if (this->GetLeftNode()!=NULL)
  {
    this->GetLeftNode()->RemoveSubTree();
    delete this->GetLeftNode();
    this->SetLeftNode(NULL);
  };

  if (this->GetRightNode()!=NULL)
  {
    this->GetRightNode()->RemoveSubTree();
    delete this->GetRightNode();
    this->SetRightNode(NULL);
  };
  return;
};

List Node::DumpDotFile()
{
  std::vector <Node *> node;
  node.push_back(this);
  int k=this->GetObservation()->GetK();
  int i;
  std::vector <std::string> var; // name of node
  std::vector <int> node_size; // size of node
  std::vector <std::vector <int> > class_table;// classification table
  std::string X="X";
  std::string N="N";
  std::vector <double> split; // split point
  std::vector <std::string> nodename;
  std::vector <std::string> leftson;
  std::vector <std::string> rightson;
  std::vector <int> index;
  index.push_back(0);
  for (i=0;i<node.size();i++)
  {
    std::string whichn = std::to_string(index[i]);
    nodename.push_back(N+whichn);
    std::vector<int> subject(node[i]->GetSubjectList());
    const Observation *obs = node[i]->GetObservation();

    int l = subject.size();
    int s;
    std::vector <int> s1;
    int j;
    for (j = 0; j < k; j++)
    {
      int count = 0;
      for (s = 0; s < l; s++)
      {
        int nr = subject[s];
        if (obs->GetY(nr)==j) count++;
      };
      s1.push_back(count);
    };

    if (std::accumulate(s1.begin(),s1.end(),0.0)!=l) {
      Rcout<<"Something wrong with node size!"<<endl;
    }
    node_size.push_back(l);
    class_table.push_back(s1);

    if (node[i]->GetRightNode()==NULL)//leaf
    {
      var.push_back("leaf");
      split.push_back(-1);
      leftson.push_back("N/A");
      rightson.push_back("N/A");

    }
    else // internal node
    {
      std::string which = std::to_string(node[i]->GetSplitVariable()+1);
      var.push_back(X+which);
      split.push_back(node[i]->GetSplitLevel());

      node.push_back(node[i]->GetLeftNode());
      index.push_back(node.size()-1);
      std::string whichleft = std::to_string(index[index.size()-1]);
      leftson.push_back(N+whichleft);

      node.push_back(node[i]->GetRightNode());
      index.push_back(node.size()-1);
      std::string whichright = std::to_string(index[index.size()-1]);
      rightson.push_back(N+whichright);

    };
  };
  return List::create(Named("node_name")=nodename,
                      Named("split_var")=var,
                      Named("node_size")=node_size,
                      Named("split")=split,
                      Named("table")=class_table,
                      Named("leftson")=leftson,
                      Named("rightson")=rightson);
};

std::vector <double> Node::GetSubjectObs()
{
  std::vector <double> Result;
  std::vector <int> sub=this->GetSubjectList();
  for (int i=0;i<sub.size();i++)
  {
    Result.push_back(this->GetObservation()->GetY(sub[i]));
  };
  return Result;
};
