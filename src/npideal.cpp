#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// #include <RcppArmadilloExtensions/sample.h>


// #include <iterator>
// #include <algorithm>
// #include <cassert>
// #include <cmath>

typedef arma::imat integerMatrixType;
typedef arma::ivec integerVectorType;
typedef arma::mat numericMatrixType;
typedef arma::vec numericVectorType;
typedef arma::uvec sizeVectorType;

template<class matrix_x>
size_t nrows(const matrix_x& iMat)
{
  return iMat.n_rows;
}

template<class matrix_x>
size_t ncols(const matrix_x& iMat)
{
  return iMat.n_cols;
}

double vectorMean(const numericVectorType& iVector)
{
  return arma::mean(iVector);
}

numericVectorType largestEigenvalue(const numericMatrixType& iMat)
{
  numericVectorType eigval;
  numericMatrixType eigvec;
  eig_sym(eigval, eigvec, iMat);
  return eigvec.col(eigvec.n_cols-1);
}


typedef std::vector<size_t> legisList;
typedef size_t int_type;

template<class int_x>
int_x unif_samp(const int_x x)
{
  return static_cast<int_x>(double(x) * unif_rand());
}

void sum_pair_votes(const integerMatrixType& v1,
		    const size_t leg1,
		    const size_t leg3, const size_t leg4,
		    size_t& v1_together, size_t& v1_apart)
{
    v1_together = 0;
    v1_apart = 0;
    for (size_t k = 0; k < nrows(v1); ++k)
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

void sum_votes(const integerMatrixType& v1,
               const size_t leg1, const size_t leg2,
               const size_t leg3, const size_t leg4,
               size_t& v1_together, size_t& v1_apart,
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

void sum_votes(const integerMatrixType& v1,
               const size_t leg1, const size_t leg2,
               const size_t leg3, const size_t leg4,
               size_t& v1_together, size_t& v1_apart)
{
    v1_together = 0;
    v1_apart = 0;
    for (size_t k = 0; k < nrows(v1); ++k)
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

// double pvalX(const int_type stat_1,
//              const int_type stat_2,
//              const int_type stat_3,
//              const int_type stat_4)
// {
//     // two-sided p-value
//     double p1 = pbinom(std::min(stat_1,stat_2),stat_1+stat_2, 0.5,  1,0);
//     // conditional minimum value for other stat
//     int_type q2a = qbinom(p1 - 1e-5,stat_3+stat_4, 0.5,  1,0);
//     q2a = std::min(q2a,stat_3+stat_4-q2a);
//     // conditional maximum value for other stat
//     int_type q2b = stat_3+stat_4-q2a;
//     // probability other less than minimum
//     double px2a = pbinom(q2a-1,stat_3+stat_4, 0.5,  1,0);
//     // probability other less than or equal to maximum
//     double px2b = pbinom(q2b,stat_3+stat_4, 0.5,  1,0);
//     // probability other less than or equal to observed
//     double px2c = pbinom((stat_1 >= stat_2) ? stat_3 : stat_4, stat_3+stat_4, 0.5,  1,0);
//     // probability other in range
//     double px2d = px2c-px2a;
//     // return probability other below given in range
//     double px2e = px2b-px2a;
//     return log(px2d/px2e);
// }


double singletest(const integerMatrixType& v1, const integerMatrixType& v2,
                  const size_t leg1, const size_t leg2,
                  const size_t leg3, const size_t leg4,
                  const legisList& voteSample1,
                  const legisList& voteSample2,
                  int_type& df)
{
  // // For testing:
  // ++df;
  // return(log(R::unif_rand()));
  
  // assert((voteSample1.size() == nrows(v1)) || (voteSample1.size() == 0));
  // assert((voteSample2.size() == nrows(v2)) || (voteSample2.size() == 0));

    size_t v1_together = 0;
    size_t v1_apart = 0;
    if (voteSample1.size() == 0)
        sum_votes(v1,leg1,leg2,leg3,leg4,v1_together,v1_apart);
    else
        sum_votes(v1,leg1,leg2,leg3,leg4,v1_together,v1_apart,voteSample1);

    size_t v2_together = 0;
    size_t v2_apart = 0;
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

    //    return Rcpp::List::create(Rcpp::Named("log.prob") = lp,
    //                              Rcpp::Named("v1.t") = v1_together,
    //                              Rcpp::Named("v1.a") = v1_apart,
    //                              Rcpp::Named("v2.t") = v2_together,
    //                              Rcpp::Named("v2.a") = v2_apart
    //                              ) ;
}

double independenttest(const integerMatrixType& v1,
                       const integerMatrixType& v2,
                       const std::vector<size_t>& legis,
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
    //    return pchisq(lp, legis.size()/4, 0, 1);
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &A, const std::vector<T> &B)
{
    std::vector<T> AB;
    AB.reserve( A.size() + B.size() );                // preallocate memory
    AB.insert( AB.end(), A.begin(), A.end() );        // add A;
    AB.insert( AB.end(), B.begin(), B.end() );        // add B;
    return AB;
}

typedef std::vector<int_type> vec_type;
std::vector<vec_type> fibyFocus(const vec_type& x, const size_t& n)
{
    std::vector<vec_type> out(0);
    if (x.size() < 4) {
        out.push_back(x);
    }
    else
      {
	size_t xs1 = x.size() - 3;
	size_t xs2 = x.size() - 2;
	size_t xs3 = x.size() - 1;
	size_t xs4 = x.size() - 0;
	if (n > 0)
	  xs1 = 1;
	if (n > 1)
	  xs2 = 2;
	if (n > 2)
	  xs3 = 3;
	if (n > 3)
	  xs4 = 4;
	for (size_t z1 = 0; z1 < xs1; ++z1)
	  for (size_t z2 = z1+1; z2 < xs2; ++z2)
            for (size_t z3 = z2+1; z3 < xs3; ++z3)
	      for (size_t z4 = z3+1; z4 < xs4; ++z4)
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

typedef std::vector<int_type> vec_type;
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
		  for (const vec_type& y1 : y)
		    out.push_back(y1+z);
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
		  for (const vec_type& y1 : y)
		    out.push_back(y1+z);
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
		  for (const vec_type& y1 : y)
		    out.push_back(y1+z);
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
		  for (const vec_type& y1 : y)
		    out.push_back(z+y1);
		  std::swap(z[1],z[2]);//1324
		  for (const vec_type& y1 : y)
		    out.push_back(z+y1);
		  std::swap(z[1],z[3]);//1423
		  for (const vec_type& y1 : y)
		    out.push_back(z+y1);
		}
      }
    return out;
}

std::vector<vec_type> fibyz(const vec_type& x)
{
    if (x.size() < 4) {
        std::vector<vec_type> y;
        y.push_back(x);
        return y;
    }
    std::vector<vec_type> out(0);
    for (size_t z1 = 0; z1 < x.size() - 3; ++z1)
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
                    for (const vec_type& y1 : y)
                        out.push_back(z+y1);
                    std::swap(z[1],z[2]);//1324
                    for (const vec_type& y1 : y)
                        out.push_back(z+y1);
                    std::swap(z[1],z[3]);//1423
                    for (const vec_type& y1 : y)
                        out.push_back(z+y1);
                  }
    return out;
}


// [[Rcpp::export]]
size_t computeCombinations(const size_t n)
{
    size_t c = 1;
    size_t i = 1;
    size_t j = 0;
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
List votetest(const integerMatrixType& v1, const integerMatrixType& v2,
	      IntegerVector focus = IntegerVector::create(),
	      size_t numLegisShuffles = 10,
	      const size_t numBootstrapSims = 10)
{
    // pcg_extras::seed_seq_from<std::random_device> seed_source;
    // rng_type rng(seed_source);

  integerVectorType focusLegis = as<integerVectorType>(focus);
  legisList legis(focusLegis.begin(), focusLegis.end());
  for (size_t i = 0; i < ncols(v1); ++i)
    if (std::find(focusLegis.begin(), focusLegis.end(), i) == focusLegis.end())
      {
      legis.push_back(i);
      }
  // for (size_t i : legis)
  //   Rcout << i << std::endl;
    // if (focusLegis.size() > 0)
    //   std::swap(legis[0],legis[focusLegis]);

    std::vector<std::vector<size_t> > legisShuffles(0);

    size_t combos = computeCombinations(legis.size());
    // Rcout << "LS: " << combos << std::endl;

    if (focusLegis.size() > 0)
      combos = 1;

    if (combos <= numLegisShuffles)
      {
	if (focusLegis.size() > 0)
	  legisShuffles = fibyFocus(legis,focusLegis.size());
	else
	  legisShuffles = fiby(legis);
        numLegisShuffles = legisShuffles.size();
        // Rcout << "DS: " << numLegisShuffles << std::endl;
        // assert(numLegisShuffles == combos);
      }
    else
      {
        for (size_t i = 0; i < numLegisShuffles; ++i)
          {
            // // std::shuffle(legis.begin(), legis.end(), rng);
            // // legisShuffles.push_back(legis);
	    // legisShuffles.push_back(RcppArmadillo::sample(legis,legis.size(),false));
	    std::random_shuffle(legis.begin() + focusLegis.size(),
				legis.end(), unif_samp<int>);
            legisShuffles.push_back(legis);
          }
        // Rcout << "RS: " << numLegisShuffles << std::endl;
      }

    legisList voteShuffle1(0);
    legisList voteShuffle2(0);

    voteShuffle1.resize(nrows(v1));
    voteShuffle2.resize(nrows(v2));
    for (size_t i = 0; i < voteShuffle1.size(); ++i)
      voteShuffle1[i] = i;
    for (size_t i = 0; i < voteShuffle2.size(); ++i)
      voteShuffle2[i] = i;

    double observedDevStat = 0.0;
    int_type df = 0, dfMax = 0;
    for (std::vector<legisList>::const_iterator ix = legisShuffles.begin();
         ix != legisShuffles.end();
         ++ix)
      {
        observedDevStat += independenttest(v1,v2,*ix,voteShuffle1,voteShuffle2,df);
        dfMax += legis.size() / 4;
      }

    // std::uniform_int_distribution<size_t> bootstrapSample1(0,nrows(v1)-1);
    // std::uniform_int_distribution<size_t> bootstrapSample2(0,nrows(v2)-1);

    voteShuffle1.resize(nrows(v1));
    voteShuffle2.resize(nrows(v2));

    double sumSimStat = 0.0;
    double sumSimStatSquared = 0.0;
    double sumSimStat2 = 0.0;
    double sumSimStat2Squared = 0.0;
    for (size_t i = 0; i < numBootstrapSims; i++)
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

    sumSimStat2 /= numBootstrapSims;
    sumSimStat2Squared /= numBootstrapSims;
    sumSimStat2Squared -= sumSimStat2*sumSimStat2;
    
    double simMean = sumSimStat / numBootstrapSims;
    double simStdDev = std::sqrt((sumSimStatSquared
				  - (sumSimStat) * (sumSimStat / numBootstrapSims))/(numBootstrapSims-1));

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
			Named("shuffles") = numLegisShuffles,
			// Named("shuf") = legisShuffles,
			Named("k") = ncols(v1),
			Named("n") = sz,
			Named("p") = pbs,
			Named("p.con") = R::pchisq(observedDevStat/numLegisShuffles, 2*df/numLegisShuffles, 0, 0));
    	
    
    // Rcout << observedDevStat << " ~ N(" << df << ", " << simStdDev << ")" << std::endl;

    // //    Rcout << pchisq(observedDevStat/numLegisShuffles, dfMax/numLegisShuffles, 0, 0) << std::endl;
    // //
    // //    Rcout << pchisq(observedDevStat, dfMax, 0, 0) << std::endl;

    // Rcout << "CM: " << R::pchisq(observedDevStat/numLegisShuffles, dfMax/numLegisShuffles, 0, 0) << std::endl;
    // //    Rcout << "BS: " << pnorm(observedDevStat,df,simStdDev,  0,0) << std::endl;
    
    // return R::pnorm(observedDevStat,df,simStdDev,  0,0);
}

template<class T>
numericVectorType rankFun(const T& x)
{
  sizeVectorType xo = sort_index(x);
  numericVectorType y(xo.size());
  size_t iLow = 0;
  for (size_t i = 0; i < xo.size(); ++i)
    {
      if ((i+1 >= xo.size()) || (x[xo[i+1]] != x[xo[i]]))
	{
	  for (size_t j = iLow; j <= i; ++j)
	    y[xo[j] ] = double(iLow+i)/2.0;
	  iLow = i+1;
	}
    }
  return y;
}

numericVectorType rankAgree(const integerMatrixType& v,
			    const size_t i,
			    const size_t j)
{
  size_t v_together, v_apart;
  numericVectorType out(ncols(v));
  size_t count = 0;
  for (size_t k = 0; k < ncols(v); ++k)
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
  numericVectorType outind = rankFun(out);
  double mean = double(ncols(v)+count+1)/2.0-1.0;
  for (size_t k = 0; k < ncols(v); ++k)
    if (out[k] < 0)
      outind[k] = mean;
  outind -= vectorMean(outind);
  return(outind);
}

// [[Rcpp::export]]
numericVectorType rankAgreeR(const integerMatrixType& v, const size_t i, const size_t j)
{
  return rankAgree(v, i-1, j-1);
}

// [[Rcpp::export]]
NumericVector voteest(const integerMatrixType& v,
		      size_t right = 0,
		      size_t left = 0,
		      const size_t maxLegisSamples = 10000)
{
  size_t n = ncols(v);
  size_t ncol = (n*(n-1))/2;
  bool useExact = ncol <= maxLegisSamples;
  if (!useExact)
    ncol = maxLegisSamples;
  numericMatrixType rmat(n,ncol);
  if (useExact)
    {
      size_t k = 0;
      for (size_t i = 1; i < ncols(v); ++i)
	for (size_t j = 0; j < i; ++j)
	  rmat.col(k++) = rankAgree(v, i, j);
    }
  else
    {
      for (size_t k = 0; k < ncol; ++k)
	{
	  size_t i = unif_samp(n-1);
	  size_t j = unif_samp(n-1);
	  if (j >= i)
	    ++j;
	  rmat.col(k) = rankAgree(v, i, j);
	}
    }
  // numericVectorType eigval;
  // numericMatrixType eigvec;
  // eig_sym(eigval, eigvec, rmat * rmat.t());
  // numericVectorType estvec = rankFun(eigvec.col(eigvec.n_cols-1))+1;
  numericVectorType estvec = rankFun(largestEigenvalue(rmat * rmat.t()))+1;
  if ((right > 0) && (right <= estvec.size()))
    {
      --right;
      if ((left > 0) && (left <= estvec.size()) && (estvec[--left] > estvec[right]))
	{
	  estvec *= -1;
	  estvec += estvec.size()+1;
	}
      else if (estvec[right]*2 < estvec.size())
       {
	  estvec *= -1;
	  estvec += estvec.size()+1;
       }
    }
    
  NumericVector est = wrap(estvec);
  est.attr("dim") = R_NilValue;
  return est;
}
