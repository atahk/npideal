#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

#define vectorMean(X)  (X.mean())
#define vectorNormalize(X)  (X.normalized())
#define nrows(X)  (X.rows())
#define ncols(X)  (X.cols())
#define matrixTrace(X)  (X.trace())
#define matrixTranspose(X) (X.transpose())
#define vectorIncrement(X,Y) X.array() += Y
#define vectorDecrement(X,Y) X.array() -= Y
typedef EIGEN_DEFAULT_DENSE_INDEX_TYPE int_type;

// For Armadillo
// #define vectorMean(X)  (arma::mean(X))
// #define vectorNormalize(X)  (arma::normalise(X))
// #define nrows(X)  (X.n_rows)
// #define ncols(X)  (X.n_cols)
// #define matrixTrace(X)  (arma::trace(X))
// #define matrixTranspose(X) (X.t())
// #define vectorIncrement(X,Y) X += Y
// #define vectorDecrement(X,Y) X -= Y
// #define largestDiag(X) (x.diag().index_max())
// typedef size_t int_type;



typedef std::vector<int_type> vec_type;
typedef vec_type legisList;


int largestDiag(const Eigen::MatrixXd& X)
{
  double max = X(0,0);
  int ind = 0;
  for (int i = 1; i < ncols(X); ++i)
    if (X(i,i) > max) {
      max = X(i,i);
      ind = i;
    }
  return ind;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &A, const std::vector<T> &B)
{
  std::vector<T> AB;
  AB.reserve( A.size() + B.size() );
  AB.insert( AB.end(), A.begin(), A.end() );
  AB.insert( AB.end(), B.begin(), B.end() );
  return AB;
}

template<class int_x>
int_x unif_samp(const int_x x)
{
  return static_cast<int_x>(double(x) * unif_rand());
}

Eigen::VectorXd largestEigenvalue(const Eigen::MatrixXd& iMat)
{
  Eigen::MatrixXd x = iMat;
  double x_trace = matrixTrace(x);
  unsigned int iter = 0;
  double threshold = 1e-8;
  do {
    x /= x_trace;
    x = matrixTranspose(x) * x;
    x_trace = matrixTrace(x);
    ++iter;
    if (iter % 10 == 0)
      threshold *= 10;
  } while (1.0 - x_trace > threshold);
  return vectorNormalize(x.col(largestDiag(x)));
}

void sum_pair_votes(const Eigen::MatrixXi& v1,
		    const int_type leg1,
		    const int_type leg3, const int_type leg4,
		    int_type& v1_together, int_type& v1_apart)
{
  v1_together = 0;
  v1_apart = 0;
  for (int_type k = 0; k < nrows(v1); ++k)
    {
      if ((v1(k,leg1)!=0) && (v1(k,leg3)*v1(k,leg4)<0))
        {
          if (v1(k,leg1)*v1(k,leg3)>0)
            ++v1_together;
          else
            ++v1_apart;
        }
    }
}

void sum_votes(const Eigen::MatrixXi& v1,
               const int_type leg1, const int_type leg2,
               const int_type leg3, const int_type leg4,
               int_type& v1_together, int_type& v1_apart,
               const legisList& voteSample)
{
  v1_together = 0;
  v1_apart = 0;
  for (legisList::const_iterator k = voteSample.begin();
       k != voteSample.end();
       ++k)
    {
      if ((v1(*k,leg1)*v1(*k,leg2)<0) && (v1(*k,leg3)*v1(*k,leg4)<0))
        {
          if (v1(*k,leg1)*v1(*k,leg3)>0)
            ++v1_together;
          else
            ++v1_apart;
        }
    }
}

void sum_votes(const Eigen::MatrixXi& v1,
               const int_type leg1, const int_type leg2,
               const int_type leg3, const int_type leg4,
               int_type& v1_together, int_type& v1_apart)
{
  v1_together = 0;
  v1_apart = 0;
  for (int_type k = 0; k < nrows(v1); ++k)
    {
      if ((v1(k,leg1)*v1(k,leg2)<0) && (v1(k,leg3)*v1(k,leg4)<0))
        {
          if (v1(k,leg1)*v1(k,leg3)>0)
            ++v1_together;
          else
            ++v1_apart;
        }
    }
}

double pval(const int_type stat_1,
            const int_type stat_2,
            const int_type stat_3,
            const int_type stat_4)
{
  double p1 = R::pbinom(std::min(stat_1,stat_2),stat_1+stat_2, 0.5,  1,1);
  int_type q2a = R::qbinom(p1 - 1e-5,stat_3+stat_4, 0.5,  1,1);
  q2a = std::min(q2a,stat_3+stat_4-q2a);
  double px2c = R::pbinom((stat_1 >= stat_2) ? stat_3 : stat_4, stat_3+stat_4, 0.5,  1,1);
  if (q2a == 0)
    return px2c;
  int_type q2b = stat_3+stat_4-q2a;
  double px2a = R::pbinom(q2a-1,stat_3+stat_4, 0.5,  1,1);
  double px2b = R::pbinom(q2b,stat_3+stat_4, 0.5,  1,1);
  double px2d = px2c+log(1.0-exp(px2a-px2c));
  double px2e = px2b+log(1.0-exp(px2a-px2b));
  return px2d-px2e;
}


double singletest(const Eigen::MatrixXi& v1, const Eigen::MatrixXi& v2,
                  const int_type leg1, const int_type leg2,
                  const int_type leg3, const int_type leg4,
                  const legisList& voteSample1,
                  const legisList& voteSample2,
                  int_type& df)
{
  int_type v1_together = 0;
  int_type v1_apart = 0;
  if (voteSample1.size() == 0)
    sum_votes(v1,leg1,leg2,leg3,leg4,v1_together,v1_apart);
  else
    sum_votes(v1,leg1,leg2,leg3,leg4,v1_together,v1_apart,voteSample1);

  int_type v2_together = 0;
  int_type v2_apart = 0;
  if (voteSample2.size() == 0)
    sum_votes(v2,leg1,leg2,leg3,leg4,v2_together,v2_apart);
  else
    sum_votes(v2,leg1,leg2,leg3,leg4,v2_together,v2_apart,voteSample2);

  double lp = 0;
  double p1, p2;

#define USE_CONDITIONAL_METHOD
#ifdef USE_CONDITIONAL_METHOD
  if ((v1_together+v1_apart>0) && (v2_together+v2_apart>0))
    {
      if ((v1_together==v1_apart) || (v2_together==v2_apart))
        {
          return(0);
        }
      if ((v1_together>v1_apart) == (v2_together>v2_apart))
        {
          return(0);
        }
      else
        {
          ++df;
          if (v1_together < v1_apart)
            std::swap(v1_together, v1_apart);
          if (v2_together < v2_apart)
            std::swap(v2_together, v2_apart);
          p1 = R::pbinom(v1_apart, v1_together+v1_apart, 0.5,  1,1)-
            R::pbinom((v1_together+v1_apart)/2,
                      v1_together+v1_apart, 0.5,  1,1);
          p2 = R::pbinom(v2_apart, v2_together+v2_apart, 0.5,  1,1)-
            R::pbinom((v2_together+v2_apart)/2,
                      v2_together+v2_apart, 0.5,  1,1);
          return(std::max(p1,p2));
        }
    }
#else

  if ((v1_together+v1_apart>0) && (v2_together+v2_apart>0))
    {
      ++df;
      p1 = R::pbinom(v2_together > v2_apart ? v1_together : v1_apart, v1_together+v1_apart, 0.5,  1,1);
      p2 = R::pbinom(v1_together > v1_apart ? v2_together : v2_apart, v2_together+v2_apart, 0.5,  1,1);
      if (v2_together == v2_apart)
        p1 = std::max(p1, R::pbinom(v1_together, v1_together+v1_apart, 0.5,  1,1));
      if (v1_together == v1_apart)
        p2 = std::max(p1, R::pbinom(v2_together, v2_together+v2_apart, 0.5,  1,1));
      return(std::max(p1,p2));
    }
#endif
  /*
    if ((v1_together+v1_apart>0) && (v2_together+v2_apart>0))
    {
    ++df;
    p1 = R::pbinom(std::min(v1_together,v1_apart), v1_together+v1_apart, 0.5,  1,1);
    p2 = R::pbinom(std::min(v2_together,v2_apart), v2_together+v2_apart, 0.5,  1,1);
    if (p1<=p2)
    lp = pval(v1_together, v1_apart, v2_together, v2_apart);
    else
    lp = pval(v2_together, v2_apart, v1_together, v1_apart);
    }
  */

  return(lp);
}

double independenttest(const Eigen::MatrixXi& v1,
                       const Eigen::MatrixXi& v2,
                       const std::vector<int_type>& legis,
                       const legisList& voteSample1,
                       const legisList& voteSample2,
                       int_type& df)
{
  double lp = 0.0;
  for (size_t i = 3; i < legis.size(); i+=4)
    lp += -2*singletest(v1, v2,
                        legis[i-3], legis[i-2], legis[i-1], legis[i],
                        voteSample1, voteSample2,
                        df);
  return lp;
}

std::vector<vec_type> fibyFocus(const vec_type& x, const int_type& n)
{
  std::vector<vec_type> out(0);
  if (x.size() < 4) {
    out.push_back(x);
  }
  else
    {
      int_type xs1 = x.size() - 3;
      int_type xs2 = x.size() - 2;
      int_type xs3 = x.size() - 1;
      int_type xs4 = x.size() - 0;
      if (n > 0)
        xs1 = 1;
      if (n > 1)
        xs2 = 2;
      if (n > 2)
        xs3 = 3;
      if (n > 3)
        xs4 = 4;
      for (int_type z1 = 0; z1 < xs1; ++z1)
        for (int_type z2 = z1+1; z2 < xs2; ++z2)
          for (int_type z3 = z2+1; z3 < xs3; ++z3)
            for (int_type z4 = z3+1; z4 < xs4; ++z4)
              {
                vec_type z(0);
                z.push_back(x[z1]);
                z.push_back(x[z2]);
                z.push_back(x[z3]);
                z.push_back(x[z4]);
                out.push_back(z);
                std::swap(z[1],z[2]);//1324
                out.push_back(z);
                std::swap(z[1],z[3]);//1423
                out.push_back(z);
              }
    }
  return out;
}

std::vector<vec_type> fiby(const vec_type& x)
{
  std::vector<vec_type> out(0);
  if (x.size() < 4) {
    out.push_back(x);
  }
  else if (x.size() % 4 == 3)
    {
      for (size_t z2 = 0; z2 < x.size() - 2; ++z2)
        for (size_t z3 = z2+1; z3 < x.size() - 1; ++z3)
          for (size_t z4 = z3+1; z4 < x.size(); ++z4)
            {
              vec_type z(0);
              z.push_back(x[z2]);
              z.push_back(x[z3]);
              z.push_back(x[z4]);
              vec_type::iterator it;
              vec_type v(x.size()+z.size());
              it=std::set_difference(x.begin(), x.end(), z.begin(), z.end(), v.begin());
              v.resize(it-v.begin());
              std::vector<vec_type> y = fiby(v);
              for (std::vector<vec_type>::const_iterator y1 = y.begin(); y1 != y.end(); ++y1)
                out.push_back((*y1)+z);
            }
    }
  else if (x.size() % 4 == 2)
    {
      for (size_t z3 = 0; z3 < x.size() - 1; ++z3)
        for (size_t z4 = z3+1; z4 < x.size(); ++z4)
          {
            vec_type z(0);
            z.push_back(x[z3]);
            z.push_back(x[z4]);
            vec_type::iterator it;
            vec_type v(x.size()+z.size());
            it=std::set_difference(x.begin(), x.end(), z.begin(), z.end(), v.begin());
            v.resize(it-v.begin());
            std::vector<vec_type> y = fiby(v);
            for (std::vector<vec_type>::const_iterator y1 = y.begin(); y1 != y.end(); ++y1)
              out.push_back((*y1)+z);
          }
    }
  else if (x.size() % 4 == 1)
    {
      for (size_t z4 = 0; z4 < x.size(); ++z4)
        {
          vec_type z(0);
          z.push_back(x[z4]);
          vec_type::iterator it;
          vec_type v(x.size()+z.size());
          it=std::set_difference(x.begin(), x.end(), z.begin(), z.end(), v.begin());
          v.resize(it-v.begin());
          std::vector<vec_type> y = fiby(v);
          for (std::vector<vec_type>::const_iterator y1 = y.begin(); y1 != y.end(); ++y1)
            out.push_back((*y1)+z);
        }
    }
  else
    {
      for (size_t z1 = 0; z1 < 1; ++z1)
        for (size_t z2 = z1+1; z2 < x.size() - 2; ++z2)
          for (size_t z3 = z2+1; z3 < x.size() - 1; ++z3)
            for (size_t z4 = z3+1; z4 < x.size(); ++z4)
              {
                vec_type z(0);
                z.push_back(x[z1]);
                z.push_back(x[z2]);
                z.push_back(x[z3]);
                z.push_back(x[z4]);
                vec_type::iterator it;
                vec_type v(x.size()+z.size());
                it=std::set_difference(x.begin(), x.end(), z.begin(), z.end(), v.begin());
                v.resize(it-v.begin());
                std::vector<vec_type> y = fiby(v);
                for (std::vector<vec_type>::const_iterator y1 = y.begin(); y1 != y.end(); ++y1)
                  out.push_back(z+(*y1));
                std::swap(z[1],z[2]);//1324
                for (std::vector<vec_type>::const_iterator y1 = y.begin(); y1 != y.end(); ++y1)
                  out.push_back(z+(*y1));
                std::swap(z[1],z[3]);//1423
                for (std::vector<vec_type>::const_iterator y1 = y.begin(); y1 != y.end(); ++y1)
                  out.push_back(z+(*y1));
              }
    }
  return out;
}

template<class T>
T computeCombinationsInternal(const T n)
{
  T c = 1;
  T i = 1;
  T j = 0;
  for (; i+3 <= n; i += 4) {
    c *= ((i * (i+1) * (i+2) * (i+3)) / (2 * 3 * 4));
    c *= 3;
    ++j;
  }
  for (; j > 0; --j)
    c /= j;
  for (j = 0; j+i <= n; ++j) {
    c *= j+i;
    c /= j+1;
  }
  return c;
}


// [[Rcpp::export]]
long int computeCombinations(const long int n)
{
  return computeCombinationsInternal<long int>(n);
}


// [[Rcpp::export]]
List votetest(const IntegerMatrix& votes1, const IntegerMatrix& votes2,
	      IntegerVector focus = IntegerVector::create(),
	      long int numLegisShuffles = 10,
	      const long int numBootstrapSims = 10)
{
  if (numLegisShuffles < 1)
    Rcpp::stop("numLegisShuffles must be positive");
  if (numBootstrapSims < 1)
    Rcpp::stop("numBootstrapSims must be positive");
  int_type iNumLegisShuffles = static_cast<int_type>(numLegisShuffles);
  int_type iNumBootstrapSims = static_cast<int_type>(numBootstrapSims);
  const Eigen::MatrixXi v1 = as<Eigen::MatrixXi>(votes1);
  const Eigen::MatrixXi v2 = as<Eigen::MatrixXi>(votes2);
  // Eigen::VectorXi focusLegis = as<Eigen::VectorXi>(focus);
  legisList focusLegis(focus.begin(), focus.end());
  legisList legis(focus.begin(), focus.end());
  for (int_type i = 0; i < ncols(v1); ++i)
    if (std::find(focusLegis.begin(), focusLegis.end(), i) == focusLegis.end())
        legis.push_back(i);
  
  // if (focusLegis.size() > 0)
  //   std::swap(legis[0],legis[focusLegis]);

  std::vector<std::vector<int_type> > legisShuffles(0);

  int_type combos = computeCombinationsInternal(legis.size());

  if (focusLegis.size() > 0)
    combos = 1;

  if (combos <= iNumLegisShuffles)
    {
      if (focusLegis.size() > 0)
        legisShuffles = fibyFocus(legis,focusLegis.size());
      else
        legisShuffles = fiby(legis);
      iNumLegisShuffles = legisShuffles.size();
      // Rcout << "DS: " << iNumLegisShuffles << std::endl;
      // assert(iNumLegisShuffles == combos);
    }
  else
    {
      for (int_type i = 0; i < iNumLegisShuffles; ++i)
        {
          // // std::shuffle(legis.begin(), legis.end(), rng);
          // // legisShuffles.push_back(legis);
          // legisShuffles.push_back(RcppArmadillo::sample(legis,legis.size(),false));
          std::random_shuffle(legis.begin() + focusLegis.size(),
                              legis.end(), unif_samp<int>);
          legisShuffles.push_back(legis);
        }
    }

  legisList voteShuffle1(0);
  legisList voteShuffle2(0);

  voteShuffle1.resize(nrows(v1));
  voteShuffle2.resize(nrows(v2));
  for (size_t i = 0; i < voteShuffle1.size(); ++i)
    voteShuffle1[i] = static_cast<int_type>(i);
  for (size_t i = 0; i < voteShuffle2.size(); ++i)
    voteShuffle2[i] = static_cast<int_type>(i);

  double observedDevStat = 0.0;
  int_type df = 0, dfMax = 0;
  for (std::vector<legisList>::const_iterator ix = legisShuffles.begin();
       ix != legisShuffles.end();
       ++ix)
    {
      observedDevStat += independenttest(v1,v2,*ix,voteShuffle1,voteShuffle2,df);
      dfMax += legis.size() / 4;
    }

  // std::uniform_int_distribution<int_type> bootstrapSample1(0,nrows(v1)-1);
  // std::uniform_int_distribution<int_type> bootstrapSample2(0,nrows(v2)-1);

  voteShuffle1.resize(nrows(v1));
  voteShuffle2.resize(nrows(v2));

  double sumSimStat = 0.0;
  double sumSimStatSquared = 0.0;
  double sumSimStat2 = 0.0;
  double sumSimStat2Squared = 0.0;
  for (int_type i = 0; i < iNumBootstrapSims; i++)
    {
      Rcpp::checkUserInterrupt();
      int_type simDf = 0;
      voteShuffle1.resize(R::rbinom(nrows(v1)+nrows(v2),
                                    double(nrows(v1))/double(nrows(v1)+nrows(v2))));
      voteShuffle2.resize(nrows(v1) + nrows(v2) - voteShuffle1.size());
      for (legisList::iterator ix = voteShuffle1.begin(); ix != voteShuffle1.end(); ++ix)
        *ix = unif_samp(nrows(v1));
      // *ix = bootstrapSample1(rng);
      for (legisList::iterator ix = voteShuffle2.begin(); ix != voteShuffle2.end(); ++ix)
        *ix = unif_samp(nrows(v2));
      // *ix = bootstrapSample2(rng);
      double simulatedDevStat = 0.0;
      for (std::vector<legisList>::const_iterator ix = legisShuffles.begin();
           ix != legisShuffles.end();
           ++ix)
        {
          simulatedDevStat += independenttest(v1,v2,*ix,voteShuffle1,voteShuffle2,simDf);
        }
      sumSimStat2 += simulatedDevStat;
      sumSimStat2Squared += simulatedDevStat * simulatedDevStat;
      simulatedDevStat -= 2*simDf;
      sumSimStat += simulatedDevStat;
      sumSimStatSquared += simulatedDevStat * simulatedDevStat;
    }

  sumSimStat2 /= iNumBootstrapSims;
  sumSimStat2Squared /= iNumBootstrapSims;
  sumSimStat2Squared -= sumSimStat2*sumSimStat2;
    
  double simMean = sumSimStat / iNumBootstrapSims;
  double simStdDev = std::sqrt((sumSimStatSquared
                                - (sumSimStat) * (sumSimStat / iNumBootstrapSims))/(iNumBootstrapSims-1));

  double pbs = 1.0;
  if (df > 0)
    pbs = R::pnorm(observedDevStat,2*df,simStdDev,  0,0);
    
  IntegerVector sz(2);
  sz[0] = nrows(v1);
  sz[1] = nrows(v2);
    
  return List::create(Named("df") = df,
                      Named("df.max") = dfMax,
                      Named("combos") = combos,
                      Named("dev") = observedDevStat,
                      Named("sim2.mean") = sumSimStat2,
                      Named("sim2.sd") = std::sqrt(sumSimStat2Squared),
                      Named("sim.mean") = simMean,
                      Named("sim.sd") = simStdDev,
                      Named("shuffles") = iNumLegisShuffles,
                      // Named("shuf") = legisShuffles,
                      Named("k") = ncols(v1),
                      Named("n") = sz,
                      Named("p") = pbs,
                      Named("p.con") = R::pchisq(observedDevStat/iNumLegisShuffles, 2*df/iNumLegisShuffles, 0, 0));
}

// typedef Eigen::VectorXd::iterator vec_it_type;
typedef double* vec_it_type;
bool rankCmp(const vec_it_type a, const vec_it_type b)
{
  return (*a) < (*b);
}

Eigen::VectorXd rankFun(Eigen::VectorXd x)
{
  std::vector<vec_it_type> p(x.size());
  std::vector<vec_it_type>::iterator pi = p.begin();
  // for (vec_it_type xi = x.begin();
  //      xi != x.end();
  //      ++xi, ++pi)
  //   *pi = xi;
  for (int j = 0; j < x.size(); ++j, ++pi)
    {
      *pi = &(x[j]);
    }
  std::sort(p.begin(), p.end(), rankCmp);
  int_type low = 0;
  int_type high = 0;
  pi = p.begin();
  std::vector<vec_it_type>::iterator pj = p.begin();
  while (pj != p.end())
    {
      while ((pj != p.end()) && (*pi == *pj))
        {
          ++pj;
          ++high;
        }
      double value = 0.5 * double(low + high - 1);
      for (std::vector<vec_it_type>::iterator pk = pi; pk != pj; ++pk)
        **pk = value;
      pi = pj;
      low = high;
    }
  return x;
}

/*
typedef arma::uvec sizeVectorType;
Eigen::VectorXd rankFun(Eigen::VectorXd x)
{
  arma::uvec xo = sort_index(x);
  Eigen::VectorXd y(xo.size());
  int_type iLow = 0;
  for (int_type i = 0; i < xo.size(); ++i)
    {
      if ((i+1 >= xo.size()) || (x[xo[i+1]] != x[xo[i]]))
	{
	  for (int_type j = iLow; j <= i; ++j)
	    y[xo[j] ] = double(iLow+i)/2.0;
	  iLow = i+1;
	}
    }
  return y;
}
*/

Eigen::VectorXd rankAgree(const Eigen::MatrixXi& v,
                          const int_type i,
                          const int_type j)
{
  int_type v_together, v_apart;
  Eigen::VectorXd out(ncols(v));
  int_type count = 0;
  for (int_type k = 0; k < ncols(v); ++k)
    {
      if ((k == i) || (k == j))
	{
	  out[k] = -1;
          ++count;
	}
      else
	{
	  sum_pair_votes(v, k, i, j, v_together, v_apart);
	  if (v_together + v_apart > 0)
	    {
	      out[k] = double(v_together)/double(v_together + v_apart);
	    }
	  else
	    {
	      out[k] = -1;
	      ++count;
	    }
	}
    }
  Eigen::VectorXd outind = rankFun(out);
  double mean = double(ncols(v)+count+1)/2.0-1.0;
  for (int_type k = 0; k < ncols(v); ++k)
    if (out[k] < 0)
      outind[k] = mean;
  vectorDecrement(outind, vectorMean(outind));
  return(outind);
}

// [[Rcpp::export]]
NumericVector rankAgreeR(const IntegerMatrix& votes, const long int i, const long int j)
{
  if ((i < 1) || (j < 1))
    Rcpp::stop("i and j are not positive");
  const Eigen::MatrixXi v = as<Eigen::MatrixXi>(votes);
  return wrap(rankAgree(v, static_cast<int_type>(i-1), static_cast<int_type>(j-1)));
}

// [[Rcpp::export]]
NumericVector voteest(const IntegerMatrix& votes,
		      long int right = 0,
		      long int left = 0,
		      const long int maxLegisSamples = 10000)
{  
  if ((right < 0) || (left < 0))
    Rcpp::stop("right and left cannot be negative");
  int_type iright = static_cast<int_type>(right);
  int_type ileft = static_cast<int_type>(left);
  int_type iMaxLegisSamples = static_cast<int_type>(maxLegisSamples);
  const Eigen::MatrixXi v = as<Eigen::MatrixXi>(votes);
  int_type n = ncols(v);
  int_type ncol = (n*(n-1))/2;
  bool useExact = ncol <= iMaxLegisSamples;
  if (!useExact)
    ncol = iMaxLegisSamples;
  Eigen::MatrixXd rmat(n,ncol);
  if (useExact)
    {
      int_type k = 0;
      for (int_type i = 1; i < ncols(v); ++i)
	for (int_type j = 0; j < i; ++j)
	  rmat.col(k++) = rankAgree(v, i, j);
    }
  else
    {
      for (int_type k = 0; k < ncol; ++k)
	{
	  int_type i = unif_samp(n-1);
	  int_type j = unif_samp(n-1);
	  if (j >= i)
	    ++j;
	  rmat.col(k) = rankAgree(v, i, j);
	}
    }
  // Eigen::VectorXd eigval;
  // Eigen::MatrixXd eigvec;
  // eig_sym(eigval, eigvec, rmat * rmat.t());
  // Eigen::VectorXd estvec = rankFun(eigvec.col(eigvec.n_cols-1))+1;
  Eigen::VectorXd estvec = rankFun(largestEigenvalue(rmat * matrixTranspose(rmat)));
  vectorIncrement(estvec, 1);
  if ((iright > 0) && (iright <= estvec.size()))
    {
      --iright;
      if ((ileft > 0) && (ileft <= estvec.size()) && (estvec[--ileft] > estvec[iright]))
	{
	  estvec *= -1;
	  vectorIncrement(estvec, estvec.size()+1);
	}
      else if (estvec[iright]*2 < estvec.size())
        {
	  estvec *= -1;
	  vectorIncrement(estvec, estvec.size()+1);
        }
    }    
  NumericVector est = wrap(estvec);
  est.attr("dim") = R_NilValue;
  return est;
}
