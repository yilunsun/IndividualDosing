#include <RcppArmadillo.h>
#include <sstream>
#include <numeric>
#include "Random.h"
#include "Model.h"
#include "Observation.h"
#include "Density.h"
#include "NodeTree.h"
using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

NodeTree::NodeTree()
{
	parent=NULL;
	this->nLeaf=0;
	this->nLevel=0;
	this->splitVariable=-1;
	left=NULL;
	right=NULL;
	return;
};

NodeTree::NodeTree(const Observation *obs)
{
	this->SetObservation(obs);
	int n=this->GetObservation()->GetN();
	this->SubjectList.resize(n);
	for (int i=0;i<n;i++)
		SubjectList[i]=i; //Initially, SubjectList has all the observations in the dataset.*/
	parent=NULL;
	left=NULL;
	right=NULL;
	splitVariable=-1;
	for (int i = 0; i < obs->GetP(); i++)
	  Splitting_var.push_back(i);
	this->nLeaf=0;
	this->nLevel=0;
	return;
};

NodeTree::~NodeTree()
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

NodeTree::NodeTree(Node &root)
{
	parent=NULL;
	this->nLeaf=0;
	this->nLevel=0;
	Node* NewNode=root.CopySubTree();
	this->SetLeftNode(NewNode->GetLeftNode());
	this->SetRightNode(NewNode->GetRightNode());
	this->SetObservation(NewNode->GetObservation());
	this->SetSplitLevel(NewNode->GetSplitLevel());
	this->SetSplitVariable(NewNode->GetSplitVariable());
	this->SetSubjectList(NewNode->GetSubjectList());
	this->SetSplittingVar(NewNode->GetSplittingVar());
	this->SetMinimumLeafSize(NewNode->GetMinimumLeafSize());
	this->SetMinLeaf(NewNode->GetMinLeaf());
	this->SetAllMinLeaf(NewNode->GetMinLeaf());
	return;
}

void NodeTree::SetParentNode(Node *node)
{
	Rcout<<"The parent for a tree can only be NULL"<<endl;
};

NodeTree* NodeTree::CopyTree(void) const
{
	NodeTree* NewTree=new NodeTree;
	NewTree->SetObservation(this->GetObservation());
	NewTree->SetSplitLevel(this->GetSplitLevel());
	NewTree->SetSplitVariable(this->GetSplitVariable());
	NewTree->SetSplittingVar(this->GetSplittingVar());
	NewTree->SetSubjectList(this->GetSubjectList());
	NewTree->SetMinimumLeafSize(this->GetMinimumLeafSize());
	NewTree->SetMinLeaf(this->GetMinLeaf());
	if (this->GetLeftNode()!=NULL || this->GetRightNode()!=NULL)
	{
		Node* pLeft=this->GetLeftNode()->CopySubTree();
		Node* pRight=this->GetRightNode()->CopySubTree();
		NewTree->SetLeftNode(pLeft);
		NewTree->SetRightNode(pRight);
	}

	return NewTree;
};

double NodeTree::Pot(const Model &modelStructure, const Model &modelVar,  const Model &modelLikelihood, Random &ran)
{
	double pot=0;
	pot =  modelStructure.Potential(*this,ran);
	pot += modelVar.Potential(*this, ran);
	pot += modelLikelihood.Potential(*this, ran);
	potential=pot;
	return pot;
};

double NodeTree::LogLikelihood(const Model &modelLikelihood, Random &ran)
{
	double pot=0;
	pot += modelLikelihood.Potential(*this, ran);
	return -pot;
};

int NodeTree::GetSize() const
{
	NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;
	int Size;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};
	Size=leaves.size();
	delete Tree;

	return Size;
};

int NodeTree::GetMiniNodeSize() const
{
  NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
  node.push_back(Tree);
  std::vector <Node *> leaves;
  
  for (int i=0;i<node.size();i++)
  {
    if (node[i]->GetRightNode()==NULL)
    {
      leaves.push_back(node[i]);
    }
    else
    {
      node.push_back(node[i]->GetLeftNode());
      node.push_back(node[i]->GetRightNode());
    };
  };
  
  int MiniNodeSize = 100000;
  
  for (int i=0;i<leaves.size();i++) //loop through leaves
  {
    std::vector <int> SubjectList=leaves[i]->GetSubjectList();
    if (MiniNodeSize >= SubjectList.size()) {
      MiniNodeSize = SubjectList.size();
    };
  };
  
  
  delete Tree;
  
  return MiniNodeSize;
};

std::vector <Node *> NodeTree::GetLeaves()
{
	NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};
	return leaves;
};

List NodeTree::GetVal() //get value of a tree and optimal label for each node
{
  double value = 0; //value is the value of the tree
  std::vector <int> label; //label is optimal dose vector for each leaf
  
  NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
  node.push_back(Tree);
  std::vector <Node *> leaves;
  int i;
  
  for (i=0;i<node.size();i++)
  {
    if (node[i]->GetRightNode()==NULL)
    {
      leaves.push_back(node[i]);
    }
    else
    {
      node.push_back(node[i]->GetLeftNode());
      node.push_back(node[i]->GetRightNode());
    };
  };
  
  for (i=0;i<leaves.size();i++) //loop through leaves
  {
    std::vector <int> SubjectList=leaves[i]->GetSubjectList();
    
    int opt = 0;
    
    // start retriving data for ith leaf
    NumericVector y_leaf(SubjectList.size()), a_leaf(SubjectList.size()), a_eval = obs->GetCanDose();
    NumericMatrix x_leaf(SubjectList.size(), obs->GetP());
    
    for (int j=0; j<SubjectList.size(); j++){ //loop within leaf node to get all observed data in that leaf
      y_leaf(j) = obs->GetV(SubjectList[j]);
      x_leaf(j,_) = obs->GetXi(SubjectList[j]);
      a_leaf(j) = obs->GetA(SubjectList[j]); // observed dose level vector
    };

    //Rcout<<"x_leaf max: "<<max(x_leaf)<<"\n";
    //Rcout<<"values: "<<value<<"\t";
    std::vector<double> est = dose_kernel_rf_cpp(y_leaf, a_leaf, x_leaf, a_eval, -999, -999, 100); //vector of estimated reward
    double maxval;
    vector_max(est, maxval, opt);
    value = value + maxval;
    if (std::isnan(value)) {
      List temp_list = List::create(Named("y_leaf")=y_leaf,
                               Named("a_leaf")=a_leaf,
                               Named("x_leaf")=x_leaf);
      return List::create(Named("label")=temp_list,
                          Named("value")=-999,
                          Named("fail")=true);
    }
      
    label.push_back(opt);
  };
  
  delete Tree;
  leaves.clear();
  node.clear();
  
  return List::create(Named("label")=label,
                      Named("value")=value,
                      Named("fail")=false);
};

std::vector <std::vector <double> > NodeTree::GetLeavesObservations()
{
  std::vector <std::vector <double> > Result;
	//vector <Node *> leaves=GetLeaves();
	//int i;

	NodeTree *Tree=this->CopyTree();
	std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};

	for (i=0;i<leaves.size();i++)
	{
	  std::vector <double> temp;
	  std::vector <int> SubjectList=leaves[i]->GetSubjectList();
		for (int j=0;j<SubjectList.size();j++)
			temp.push_back(obs->GetY(SubjectList[j]));
		Result.push_back( temp );
	};
	delete Tree;
	leaves.clear();
	node.clear();


	return Result;
};

std::vector <std::vector <int> > NodeTree::GetLeavesSubjects()
{
  std::vector <std::vector <int> > Result;
	//vector <Node *> leaves=GetLeaves();
	//int i;

	NodeTree *Tree=this->CopyTree();
	std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{

			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};

	for (i=0;i<leaves.size();i++)
	{
	  std::vector <int> SubjectList=leaves[i]->GetSubjectList();
		Result.push_back(SubjectList);
	};
	delete Tree;
	leaves.clear();
	node.clear();

	return Result;
};

std::vector <int> NodeTree::GetSplittingVariableSet() const
{
  std::vector <int> Result;
	//vector <Node *> leaves=GetLeaves();
	//int i;

	NodeTree *Tree=this->CopyTree();
	std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;

	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL)
		{
			leaves.push_back(node[i]);
		}
		else
		{
			Result.push_back(node[i]->GetSplitVariable());
			node.push_back(node[i]->GetLeftNode());
			node.push_back(node[i]->GetRightNode());
		};
	};

	delete Tree;
	leaves.clear();
	node.clear();

	return Result;
};

bool NodeTree::Misclassfication(int &LeavesSize, int &Error, Density *like) const
{

	NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
	node.push_back(Tree);
	std::vector <Node *> leaves;
	int i;
	int k = Tree->GetObservation()->GetK();

	LeavesSize=0; Error=0;


	for (i=0;i<node.size();i++)
	{
		if (node[i]->GetRightNode()==NULL) //leaf
		{
			leaves.push_back(node[i]);

		  std::vector<int> subject(node[i]->GetSubjectList());
			const Observation *obs = node[i]->GetObservation();
			std::vector <int> s1;

			int l = subject.size();

			int s,j;

			NumericVector leaf;

			for (s = 0; s < l; s++)
			{
			  int nr = subject[s];
			  leaf.push_back(obs->GetY(nr));
			  for (j = 0; j < k; j++){
			    int count = 0;
			    if (obs->GetY(nr)==j) count++;
			    s1.push_back(count);
			  };
			};

			std::vector <double> pMean=like->PosteriorMean(leaf);

			pMean.push_back(1-std::accumulate(pMean.begin(), pMean.end(), 0.0));
			int pmax=std::distance(pMean.begin(), std::max_element(pMean.begin(), pMean.end()));
			if (pmax==k)
			  s1.erase(s1.begin());
			else
			  s1.erase(s1.begin()+pmax-1);

			Error+=std::accumulate(s1.begin(), s1.end(), 0);
		}
		else
		{
			node.push_back(node[i]->GetRightNode());
			node.push_back(node[i]->GetLeftNode());
		};
	};
	LeavesSize=leaves.size();
	delete Tree;
	Tree=NULL;
	return true;
};

arma::colvec NodeTree::repc(const NumericVector& input, const int& ct){
  
  arma::colvec input_a = as<arma::colvec>(input);
  
  arma::colvec result = input_a;
  
  for (int i=0;i<ct-1;i++){
    result = arma::join_cols(result, input_a);
  }
  
  return result;
}

NumericVector NodeTree::seq_length(const double& start, const double& end, const int& length_out){
  NumericVector result;
  double increment = (end-start)/(length_out - 1);
  for (int i = 0; i < length_out; i++){
    result.push_back(start+i*increment);
  }
  return result;
}

NumericVector NodeTree::dirtyGBM(const NumericMatrix& x, const NumericVector& y, const NumericMatrix& newdf){ // this implementation requires expose predict in xgboost package as predict2
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:xgboost2"); 
  
  // Make function callable from C++
  Rcpp::Function xgbm = package_env["xgboost2"];  
  Rcpp::Function predict2 = package_env["predict2"]; 
  
  // Call the function and receive output (might not be list)
  //Rcout << "here";
  SEXP rf_obj = xgbm(Named("data")=x, _["label"]=y, _["max.depth"] =3 , _["eta"] = 1, _["nthread"] = 2, _["nrounds"] = 2, _["verbose"] = 0);
  //Rcout << "2";
  NumericVector yhat = predict2(rf_obj, newdf);
  
  return yhat;
}

List NodeTree::dirtyDensity(const arma::colvec& x){ // this implementation requires expose predict in xgboost package as predict2
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:stats"); 
  
  // Make function callable from C++
  Rcpp::Function densi = package_env["density"];
  
  //Rcout << "here";
  List dens = densi(x);
  
  return dens;
}

arma::colvec NodeTree::dirtySpline(const arma::colvec& x, const arma::colvec& y, const arma::colvec& newdf){ // this implementation requires expose predict in xgboost package as predict2
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:stats"); 
  
  // Make function callable from C++
  Rcpp::Function smoothspline = package_env["smooth.spline"];  
  Rcpp::Function predict2 = package_env["predict"]; 
  
  // Call the function and receive output (might not be list)
  //Rcout << "here";
  SEXP rf_obj = smoothspline(Named("x")=x, _["y"]=y);
  //Rcout << "2";
  List yhat = predict2(rf_obj, newdf);
  
  return yhat["y"];
}

arma::colvec NodeTree::w_fn(const double& bw, const arma::colvec a_vals, const arma::colvec& a, const int& n) {
  arma::colvec w_avals(a_vals.n_elem);
  for (int i=0; i<a_vals.n_elem; i++) {
    arma::colvec a_std = (a - a_vals(i))/bw;
    arma::colvec kern_std = arma::normpdf(a_std)/bw;
    w_avals(i) = arma::mean(pow(a_std,2.0) % kern_std) * (arma::normpdf(0.0)/bw)/(arma::mean(kern_std) * arma::mean(pow(a_std,2.0) % kern_std) - pow(arma::mean(a_std % kern_std),2.0));
  }
  return w_avals/n;
}

arma::colvec NodeTree::dirtyLocpoly(const double& h, const arma::colvec& out, const arma::colvec& a, const arma::colvec& a_out){
  
  Rcpp::Environment package_env("package:KernSmooth"); 
  
  Rcpp::Function locpoly = package_env["locpoly"];
  // 
  List polybw = locpoly(Named("x")=a, _["y"]=out, _["bandwidth"] = h);
  
  arma::colvec result;
  
  arma::interp1(as<arma::colvec>(polybw["x"]), as<arma::colvec>(polybw["y"]), a_out, result, "linear");
  
  return result;
}

double NodeTree::risk_fn(const double& h, const arma::colvec& a_vals, const arma::colvec& out, const arma::colvec& a, const int& n) {
  
  arma::colvec hats;
  
  arma::colvec a_interp = a;
  a_interp.elem( find(a_interp > max(a_vals)) ).fill(max(a_vals));
  a_interp.elem( find(a_interp < min(a_vals)) ).fill(min(a_vals));
  
  arma::interp1(a_vals, w_fn(h, a_vals, a, n), a_interp, hats, "linear");
  
  arma::colvec cts_eff_fn = dirtyLocpoly(h, out, a, a);
  // 
  arma::colvec one1(hats.n_elem);
  one1.ones();
  
  return arma::mean(pow((out - cts_eff_fn)/(one1 - hats),2.0));
}

std::vector<double> NodeTree::dose_kernel_rf_cpp(const NumericVector& y, const NumericVector& a, const NumericMatrix& x, const NumericVector& a_out,  
                                       double a_min = -999, double a_max = -999, const int n_pts = 100) 
{
  
  //require(KernSmooth)
  
  //kern <- function(t) {
  //  dnorm(t)
  //}
  int n = x.nrow();
  
  NumericVector bw_seq = seq_length(0.2, 10, 100);
  
  if (a_min == -999) 
  {
    a_min = min(a);
  }
  if (a_max == -999) 
  {
    a_max = max(a);
  }
  
  NumericVector a_vals = seq_length(a_min, a_max, n_pts);
  
  NumericVector pimod_vals = dirtyGBM(x, a, x);
  
  NumericVector pi2mod_vals = dirtyGBM(x, Rcpp::pow(a - pimod_vals, 2.0), x);
  
  arma::mat x_p = as<arma::mat>(x);
  
  arma::mat x_new = x_p;
  
  for (int i=0; i<n_pts; i++){
    x_new = arma::join_cols(x_new, x_p);
  }
  
  arma::colvec a_new(a_vals.size()*n);
  for (int i=0;i<a_vals.size();i++){
    for (int j=0;j<n;j++){
      a_new(i*n+j) = a_vals(i);
    }
  }
  
  arma::mat xa_new = arma::join_rows(x_new, arma::join_cols(as<arma::colvec>(a), a_new));
  
  arma::colvec muhat_vals = as<arma::colvec>(dirtyGBM(Rcpp::cbind(x,a), y, wrap(xa_new)));
  //Rcout<<"here aaa1 \n";
  arma::colvec a_std = (xa_new.col(xa_new.n_cols-1) - repc(pimod_vals,n_pts+1))/sqrt(abs(repc(pi2mod_vals,n_pts+1)));
  //Rcout<<"here aaa2 \n";
  arma::colvec pihat_vals;
  
  List density = dirtyDensity(a_std);
  
  arma::interp1(as<arma::colvec>(density["x"]), as<arma::colvec>(density["y"]), a_std, pihat_vals, "linear");
  //Rcout<<"here aaa3 \n";
  arma::colvec pihat = pihat_vals.subvec(0,n-1);
  
  arma::mat pihat_mat = arma::reshape(pihat_vals.subvec(n, pihat_vals.n_elem-1), n, n_pts);
  
  arma::colvec varpihat = dirtySpline(as<arma::colvec>(a_vals), arma::conv_to<arma::colvec>::from(arma::mean(pihat_mat, 0)), as<arma::colvec>(a));
  //Rcout<<"here aaa4 \n";
  arma::mat varpihat_mat;
  
  for (int i=0;i<n;i++){
    varpihat_mat.insert_rows(i,arma::conv_to<arma::colvec>::from(arma::mean(pihat_mat, 0)));
  }
  
  arma::colvec muhat = muhat_vals.subvec(0,n-1);
  
  arma::mat muhat_mat = arma::reshape(muhat_vals.subvec(n, muhat_vals.n_elem-1), n, n_pts);
  
  arma::colvec mhat = dirtySpline(as<arma::colvec>(a_vals), arma::conv_to<arma::colvec>::from(arma::mean(muhat_mat, 0)), as<arma::colvec>(a));
  //Rcout<<"here aaa5 \n";
  arma::mat mhat_mat;
  
  for (int i=0;i<n;i++){
    mhat_mat.insert_rows(i, arma::conv_to<arma::colvec>::from(arma::mean(muhat_mat, 0)));
  }
  
  arma::colvec pseudo_out = (as<arma::colvec>(y) - muhat)/(pihat/varpihat) + mhat;
  //Rcout<<"here aaa6 \n";
  arma::colvec risk_est;
  risk_est.ones(bw_seq.size());
  //Rcout<<"size of bw_seq: "<<bw_seq.size()<<"\n";
  //Rcout<<"size of risk_est: "<<risk_est.size()<<"\n";
  for (int i=0;i<risk_est.size();i++){
    risk_est(i) = risk_fn(bw_seq(i), a_vals, pseudo_out, as<arma::colvec>(a), n);
  }  
  //Rcout<<"here aaa7 \n";
  double h_opt = bw_seq(index_min(risk_est));
  
  arma::colvec est = dirtyLocpoly(h_opt, pseudo_out, a, a_out);
  
  return arma::conv_to<std::vector<double> >::from(est);
}


void NodeTree::vector_max(std::vector<double> v, double &max, int &imax){
  std::vector<double>::size_type p=0;
  imax = -1;
  max = std::numeric_limits<double>::lowest();
  
  for (auto &val : v)
  {
    if (!std::isnan(val) && val>max)
    {
      imax = p;
      max = val;
    }
    p++;
  }
}

void NodeTree::checkTree(void) const
{
  NodeTree *Tree=this->CopyTree();
  std::vector <Node *> node;
  node.push_back(Tree);
  std::vector <Node *> leaves;
  
  for (int i=0;i<node.size();i++)
  {
    Rcout<<"Now Checking MinLeaf is: "<<node[i]->GetMinLeaf()<<"\t";
    if (node[i]->GetRightNode()==NULL)
    {
      leaves.push_back(node[i]);
    }
    else
    {
      node.push_back(node[i]->GetLeftNode());
      node.push_back(node[i]->GetRightNode());
    };
  };
  Rcout<<"Finished Checking MinLeaf\n";
  delete Tree;
};

