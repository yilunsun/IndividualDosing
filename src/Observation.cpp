#include "Random.h"
#include "Observation.h"

//constructor
Observation::Observation(const NumericMatrix &x, const IntegerVector &y, const int &cat_num)
{
	this->x = x;
  this->y = y;
  this->cat_num = cat_num;

  this->n = x.nrow();
  this->p = x.ncol();

  this->k = max(y) + 1;
};


void Observation::Show()
{
	Rcout<<"Y:\n";
	for (int i=0; i<y.size(); ++i)
		Rcout<<y[i]<<"\t";

	Rcout<<"\n X:\n";
	for (int i=0;i<n;i++)
	{
		Rcout<<i<<":\t";
		for (int j=0;j<p;j++)
		{
			Rcout<<x(i,j)<<"\t";
			if (x(i,j)>1.0) Rcout<<"ERROR";
		};
	};
	return;
};

void Observation::SummaryStat()
{
  Rcout<<"Dimension of input data X:\n"<<n<<" by "<<p<<"\n";

  Rcout<<"Number of categories in outcome:\n"<<k<<"\n";

  Rcout<<"Number of categorical variables in X:\n"<<cat_num<<"\n";

  return;
};


Observation::~Observation()
{
  return;
};
