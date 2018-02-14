


#include "Random.h"


void Random::setseed()
{
	srand((unsigned) time(NULL));
};

Random::Random(unsigned int seed)
{
  seedValue = seed;

  haveNorm01 = 0;

  return;
};

Random::Random()
{
	setseed();
};

double Random::Unif()
{
	  return 1.0*rand()/RAND_MAX;
};


Random::~Random(void)
{
  return;
};



double Random::Unif01(void)
{
  seedValue = MULTIPLIER * seedValue + SHIFT;
  double r = ((double) seedValue) * INVMOD;

  return r;
};

double Random::PotentialUnif(double x, double lower, double upper)
{
  if (x<lower || x>upper) Rcout<<"Random::Unif is wrong!";

  double pot = - log(1/(upper-lower));

  return pot;
};

double Random::PotentialUnifDiscrete(double x, int lower, int upper)
{
  if (x<lower || x>upper) Rcout<<"Random::UnifDiscrete is wrong!";

  double pot = - log(1/(upper-lower+1));

  return pot;
};


double Random::Norm01(void)
{
  if (haveNorm01 == 1)
    {
      haveNorm01 = 0;
      return norm;
    }
  else
    {
      double u1,u2;
      double x,y;

      u1 = Unif01();
      u2 = Unif01();

      x = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
      y = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

      haveNorm01 = 1;
      norm = x;

      return y;
    }
};




double Random::Exponential(double lambda)
{
  double u = Unif01();

  while(u==0)
      u=Unif01();

  double v = - log(u) / lambda;

  return v;
};


double Random::Gamma(int alpha, double beta)
{
	// 1/gamma(alpha) * beta^alpha x^{alpha-1} exp(-x*beta);
	// = 1/alpha \sum exp(beta);
	double sum=0;
	for (int i=0;i<alpha;i++)
		sum+=this->Exponential(beta);
	return sum;
};

int Random::Poisson(double lambda)
{
  //
  // Count the number of events before 1 in a
  // Poisson process with intensity lambda.
  // Note: this is very inefficient if lambda is large!!!
  //

  int v = 0;
  double sum = Exponential(lambda);
  while (sum <= 1.0)
    {
      v++;
      sum += Exponential(lambda);
    }

  return v;
};




int Random::Binomial(int n,double p)
{
  //
  // Count the number of successes.
  // Note: this is very inefficient if n is large
  //

  int v = 0;
  int k;
  for (k = 0; k < n; k++)
    v += (Unif01() <= p);

  return v;
};



int Random::Discrete(const std::vector<double> &prob)
{
  double sum = 0.0;
  int i;
  for (i = 0; i < prob.size(); i++)
    sum += prob[i];

  double u = sum * Unif01();

  int v;
  sum = prob[0];
  for (v = 0; u > sum; v++) sum += prob[v+1];

  return v;
};




double Random::PotentialGaussian(double variance,double mean,double x)
{
  double pot = 0.5 * (x - mean) * (x - mean) / variance;
  pot += 0.5 * log(2.0 * M_PI);
  pot += 0.5 * log(variance);

  return pot;
};

double Random::PotentialExp(double lambda, double x)
{
	double pot = - (log(lambda) - lambda*x);
	return pot;
};

double Random::PotentialWeibull(double alpha,double beta, double x)
{
	double pot=0;
	pot = -(log(alpha) + log(beta) + (alpha-1)*log(x) - exp(alpha*log(x))*beta);
	return pot;
};


double Random::PotentialTruncatedGaussian(double variance,double mean,double x, double lower, double upper)
{
  if (x<lower || x>upper) Rcout<<"Random::TruncatedGaussian is wrong!";
  double c=normcdf(upper, mean, sqrt((double) variance))- normcdf(lower, mean, sqrt((double) variance));

  double pot = 0.5 * (x - mean) * (x - mean) / variance;
  pot += 0.5 * log(2.0 * M_PI);
  pot += 0.5 * log(variance);
  pot += -log(c);
  return pot;
};


double Random::PotentialPoisson(double lambda,int x)
{
  double pot;

  pot = - x * log(lambda);
  pot += lambda;
  pot += lnGamma(x + 1.0);

  return pot;
};

double Random::PotentialTruncatedPoisson(double lambda,int x, int length)
{
  std::vector <double> Probs(length);

	int i;
	double sum=0;

	for (i=0;i<length;i++)
	  //		Probs[i]=exp(-PotentialPoisson(lambda, i+1));
	  Probs[i]=exp(-PotentialPoisson(lambda, i));

	for (i=0;i<length;i++)
		sum+=Probs[i];
	for (i=0;i<length;i++)
		sum=Probs[i]/sum;
	return -log(Probs[x-1]);
};



double Random::PotentialBinomial(int n,double p,int x)
{
  double pot = - x * log(p) - (n - x) * log(1.0 - p);
  pot += - lnGamma(n + 1.0) + lnGamma(x + 1.0) + lnGamma(n - x + 1.0);

  return pot;
};



double Random::lnGamma(double xx)
{
  //
  // return ln(Gamma(x)). Numerical routine from Numerical Recipies.
  //

  double x,y,tmp,ser;
  static double cof[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
			  -0.5395239384953e-5};
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++) ser += cof[j]/++y;

  return - tmp + log(2.5066282746310005 * ser / x);
};












///////////////////////////////////////////////////////////////////////
// Some private functions
//////////////////////////////////////////////////////////////////////


double Random::erfcore(double x)
{
  std::vector <double> a(5), b(4), c(9), d(8), p(6), q(5);
  double result=0;
  a[0]=3.16112374387056560e00; a[1]=1.13864154151050156e02;
  a[2]=3.77485237685302021e02; a[3]=3.20937758913846947e03;
  a[4]=1.85777706184603153e-1;

  b[0]=2.36012909523441209e01; b[1]=2.44024637934444173e02;
  b[2]=1.28261652607737228e03; b[3]=2.84423683343917062e03;

  c[0]=5.64188496988670089e-1; c[1]=8.88314979438837594e00;
  c[2]=6.61191906371416295e01; c[3]=2.98635138197400131e02;
  c[4]=8.81952221241769090e02; c[5]=1.71204761263407058e03;
  c[6]=2.05107837782607147e03; c[7]=1.23033935479799725e03;
  c[8]=2.15311535474403846e-8;

  d[0]=1.57449261107098347e01; d[1]=1.17693950891312499e02;
  d[2]=5.37181101862009858e02; d[3]=1.62138957456669019e03;
  d[4]=3.29079923573345963e03; d[5]=4.36261909014324716e03;
  d[6]=3.43936767414372164e03; d[7]=1.23033935480374942e03;

  p[0]=3.05326634961232344e-1; p[1]=3.60344899949804439e-1;
  p[2]=1.25781726111229246e-1; p[3]=1.60837851487422766e-2;
  p[4]=6.58749161529837803e-4; p[5]=1.63153871373020978e-2;

  q[0]=2.56852019228982242e00; q[1]=1.87295284992346047e00;
  q[2]=5.27905102951428412e-1; q[3]=6.05183413124413191e-2;
  q[4]=2.33520497626869185e-3;

  if (abswu((double) x)<=0.46875)
    {
      double y=abswu((double) x);
      double z=y*y;
      double xnum=a[4]*z;
      double xden=z;
      for (int i=0;i<3;i++)
	{
	  xnum=(xnum+a[i])*z;
	  xden=(xden+b[i])*z;
	}
      result=x*(xnum+a[3])/(xden+b[3]);
    };
  if (abswu((double) x)>0.46875 && abswu((double) x)<=4)
    {
      double y=abswu((double) x);
      double xnum=c[8]*y;
      double xden=y;
      double z, del;
      for (int i=0;i<7;i++)
	{
	  xnum=(xnum+c[i])*y;
	  xden=(xden+d[i])*y;
	};
      result=(xnum+c[7])/(xden+d[7]);
      z=fix(y*16)/16;
      del=(y-z)*(y+z);
      result=exp(-z*z)*exp(-del)*result;
    };
  if (abswu((double) x)>4)
    {
      double y=abswu((double) x);
      double z=1/(y*y);
      double xnum=p[5]*z;
      double xden=z;
      double del;
      for (int i=0;i<4;i++)
	{
	  xnum=(xnum+p[i])*z;
	  xden=(xden+q[i])*z;
	};
      result=z*(xnum+p[4])/(xden+q[4]);
      result=(1/sqrt(M_PI)-result)/y;
      z=fix(y*16)/16;
      del=(y-z)*(y+z);
      result=exp(-z*z)*exp(-del)*result;
    };
  if (x>0.46875)
    result=(0.5-result)+0.5;
  if (x<-0.46875)
    result=(-0.5+result)-0.5;
  return result;
};

double Random::fix(double x)
{

  if (x>0) return floor(x);
  if (x<0) return floor(x)+1;
  return 0;
};

double Random::normcdf(double x, double mu, double sd)
{
  x=((x-mu)/sd);
  return 0.5*(1-erfcore(-x/sqrt(2.0)));
};


double Random::erf(double x)
{
  return erfcore(x);
};

double Random::erfinv(double y)
{
  double x=0, y0=0.7;
  double u=0;
  std::vector <double> a(4), b(4), c(4), d(2);

  a[0]=0.886226899;  a[1]=-1.645349621; a[2]=0.914624893;  a[3]=-0.140543331;
  b[0]=-2.118377725; b[1]=1.442710462;  b[2]=-0.329097515; b[3]=0.012229801;
  c[0]=-1.970840454; c[1]=-1.624906493; c[2]=3.429567803;  c[3]=1.641345311;
  d[0]=3.543889200;  d[1]=1.637067800;

  if (abswu((double) y)<=y0)
    {
      double z=y*y;
      x=y*( ( ( a[3]*z+a[2])*z+a[1])*z+a[0])/((((b[3]*z+b[2])*z+b[1])*z+b[0])*z+1);
    };

  if (y>y0 && y<1)
    {
      double z=sqrt(-log((1-y)/2));
      x=(((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1);
    };
  if (y<-y0 && y>-1)
    {
      double z=sqrt(-log((1+y)/2));
      x=-(((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1);
    };

  u=(erf(x)-y)/(2/sqrt(M_PI)*exp(-x*x));
  x=x-u/(1+x*u);
  return x;
};


double Random::abswu(double x)
{
	if (x<0) x=-x;
	return x;
};





double Random::norminv(double p, double mu, double sd)
{
  double x=-sqrt(2.0)*erfinv(1-2*p);
  return x*sd+mu;
};

