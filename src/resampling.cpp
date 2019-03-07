#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;


// function to permute the ancestry vector so that A^i = i for as many i as possible
void permute(IntegerVector & a, const int & nparticles){
  int swap;
  for (int i = 0; i < nparticles; i++){
    if (a(i) != i && a(a(i)) != a(i)){
      swap = a(a(i));
      a(a(i)) = a(i);
      a(i) = swap;
      i--;
    }
  }
}

// [[Rcpp::export]]
IntegerVector systematic_resampling_(int nsamples, const NumericVector & weights){
  RNGScope scope;
  IntegerVector ancestors(nsamples);
  NumericVector u_vec = runif(1);
  double u = u_vec(0) / nsamples;
  int j = 0;
  double csw = weights(0);
  for(int k = 0; k < nsamples; k++){
    while(csw < u){
      j++;
      csw += weights(j);
    }
    ancestors(k) = j;
    u = u + 1./nsamples;
  }
  return ancestors + 1; // add 1 so that smallest is 1, and largest is N
}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
IntegerVector multinomial_resampling_(int nsamples, const NumericVector & weights){
  RNGScope scope;
  IntegerVector ancestors(nsamples);
  int nparents = weights.size();
  NumericVector cumsumw = cumsum(weights);
  NumericVector uniforms = runif(nsamples);
  double sumw = cumsumw(nparents - 1);
  double lnMax = 0;
  int j = nparents;
  for (int i = nsamples; i > 0; i--){
    if (j > 0){
      lnMax += log(uniforms(i-1)) / i;
      uniforms(i-1) = sumw * exp(lnMax); // sorted uniforms
      while (j > 0 && uniforms(i-1) < cumsumw(j-1)){
        j --;
      }
      ancestors(i-1) = j;
    } else {
      ancestors(i-1) = 0;
    }
  }
  std::random_shuffle(ancestors.begin(), ancestors.end(), randWrapper);
  return ancestors + 1;
}

// [[Rcpp::export]]
IntegerVector SSP_resampling_(int nsamples, const NumericVector & weights){
  RNGScope scope;
  int N = weights.size();
  int n,m,k;
  NumericVector uniforms = runif(N); // vector of N unif(0,1)
  const double tol = 1e-15; // tolerance used to test whether a float is an integer
  NumericVector Y(N);
  IntegerVector nb_offsprings(N);
  IntegerVector ancestors(nsamples);
  for(int i = 0; i < N; i++){
    Y(i) = nsamples*weights(i);
    nb_offsprings(i) = floor(Y(i));
  }
  n = 0;
  m = 1;
  for(k = 0; k < N; k++){
    if ((m+1) > N){
      break;
    }
    double delta_n = nb_offsprings(n)+1-Y(n);
    double delta_m = Y(m)-nb_offsprings(m);
    double epsilon_n = Y(n)-nb_offsprings(n);
    double epsilon_m = nb_offsprings(m)+1-Y(m);
    //// Test if Y(n) or Y(m) is an integer
    bool Yn_is_integer = ((delta_n < tol) || (epsilon_n < tol));
    bool Ym_is_integer = ((delta_m < tol) || (epsilon_m < tol));
    //// Deal with easy cases first: either Y(n) or Y(m) (or both) is already an integer
    if (Yn_is_integer && Ym_is_integer){
      n = m+1;
      m = m+2;
    } else if (Yn_is_integer){
      n = m;
      m = m+1;
    } else if (Ym_is_integer){
      m = m+1;
    } else {
      double delta = min(delta_n, delta_m);
      double epsilon = min(epsilon_n, epsilon_m);
      //// Keeping track of delta_n, delta_m, epsilon_n, epsilon_m prevents rounding errors
      if (uniforms(k) < epsilon/(epsilon+delta)){
        if (std::abs(delta_n-delta_m)<tol){
          // Test delta_n == delta_m while accounting for numerical rounding imprecisions
          nb_offsprings(n) += 1;
          n = m+1;
          m = m+2;
        } else if (delta_n < delta_m){
          Y(m) = Y(m) - delta;
          nb_offsprings(n) += 1;
          n = m;
          m = m+1;
        } else if (delta_n > delta_m){
          Y(n) = Y(n) + delta;
          m = m+1;
          // Note that nb_offsprings(m) is already set at the correct value
        }
      } else {
        if (std::abs(epsilon_n-epsilon_m)<tol){
          // Test epsilon_n == epsilon_m while accounting for numerical rounding imprecisions
          nb_offsprings(m) += 1;
          n = m+1;
          m = m+2;
        } else if (epsilon_n < epsilon_m){
          Y(m) = Y(m) + epsilon;
          n = m;
          m = m+1;
          // Note that nb_offsprings(n) is already set at the correct value
        } else if (epsilon_n > epsilon_m){
          Y(n) = Y(n) - epsilon;
          nb_offsprings(m) += 1;
          m = m+1;
        }
      }
    }
  }
  int idx = 0;
  for (int l = 0; l < N; l++){
    for (int j = 0; j < nb_offsprings(l); j++){
      ancestors(idx) = l;
      idx++;
    }
  }
  return (ancestors + 1); // shift the indices from C++ to R convention
}
