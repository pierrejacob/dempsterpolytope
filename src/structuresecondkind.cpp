#include <RcppEigen.h>
using namespace Rcpp;

// the following draws A_n uniformly in pi_k(theta), for all n
// and also computes the minimum of A_n,ell / A_n,k  for all n,ell

// [[Rcpp::export]]
List runif_piktheta_cpp(int n, int k, NumericVector theta){
  RNGScope rngScope;
  int K = theta.length();
  NumericVector w = rexp(n * K, 1.0);
  NumericMatrix A(n, K);
  double sumrow = 0.;
  double temp = 0.;
  NumericVector minratios(K);
  std::fill(minratios.begin(), minratios.end(), 0);
  for (int irow = 0; irow < n; irow ++){
    sumrow = 0.;
    for (int icol = 0; icol < K; icol ++){
      A(irow,icol) = w(irow*K + icol);
      sumrow += A(irow,icol);
    }
    for (int icol = 0; icol < K; icol ++){
      A(irow,icol) /= sumrow;
      // at this point A(irow,icol) is Dirichlet (1,..,1)
    }
    for (int icol = 0; icol < K; icol ++){
      if (icol+1 != k){
        A(irow,icol) = A(irow,k-1) * theta(icol) + A(irow,icol);
      }
    }
    A(irow,k-1) *= theta(k-1);
    // at this point we have A_n
    // update minimum of A_n,ell/ A_n,k
    if (irow == 0){
      for (int icol = 0; icol < K; icol ++){
        minratios(icol) = A(0,icol) / A(0,k-1);
      }
    } else {
      for (int icol = 0; icol < K; icol ++){
        temp = A(irow,icol) / A(irow,k-1);
        if (temp < minratios(icol)){
          minratios(icol) = temp;
        }
      }
    }
  }
  return List::create(Named("pts") = A, Named("minratios") = minratios);
}