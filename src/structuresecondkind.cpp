#include <RcppEigen.h>
using namespace Rcpp;

// the following draws A_n uniformly in pi_k(theta), for all n
// and also computes the minimum of A_n,ell / A_n,k  for all n,ell

// [[Rcpp::export]]
List runif_piktheta_cpp(int n, int k, NumericVector theta){
  RNGScope rngScope;
  int K = theta.length();
  GetRNGstate();
  NumericVector w = rexp(n * K, 1.0);
  PutRNGstate();
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

// Coupled version, with common random numbers
// [[Rcpp::export]]
List crng_runif_piktheta_cpp(int n, int k, NumericVector & theta1, NumericVector & theta2){
  RNGScope rngScope;
  int K = theta1.length();
  GetRNGstate();
  NumericVector w = rexp(n * K, 1.0);
  PutRNGstate();
  NumericMatrix A1(n, K);
  NumericMatrix A2(n, K);
  double sumrow1 = 0.;
  double temp1 = 0.;
  NumericVector minratios1(K);
  std::fill(minratios1.begin(), minratios1.end(), 0);
  double sumrow2 = 0.;
  double temp2 = 0.;
  NumericVector minratios2(K);
  std::fill(minratios2.begin(), minratios2.end(), 0);
  for (int irow = 0; irow < n; irow ++){
    sumrow1 = 0.;
    sumrow2 = 0.;
    for (int icol = 0; icol < K; icol ++){
      A1(irow,icol) = w(irow*K + icol);
      sumrow1 += A1(irow,icol);
      A2(irow,icol) = w(irow*K + icol);
      sumrow2 += A2(irow,icol);
    }
    for (int icol = 0; icol < K; icol ++){
      A1(irow,icol) /= sumrow1;
      A2(irow,icol) /= sumrow2;
      // at this point A(irow,icol) is Dirichlet (1,..,1)
    }
    for (int icol = 0; icol < K; icol ++){
      if (icol+1 != k){
        A1(irow,icol) = A1(irow,k-1) * theta1(icol) + A1(irow,icol);
        A2(irow,icol) = A2(irow,k-1) * theta2(icol) + A2(irow,icol);
      }
    }
    A1(irow,k-1) *= theta1(k-1);
    A2(irow,k-1) *= theta2(k-1);
    // at this point we have A_n
    // update minimum of A_n,ell/ A_n,k
    if (irow == 0){
      for (int icol = 0; icol < K; icol ++){
        minratios1(icol) = A1(0,icol) / A1(0,k-1);
        minratios2(icol) = A2(0,icol) / A2(0,k-1);
      }
    } else {
      for (int icol = 0; icol < K; icol ++){
        temp1 = A1(irow,icol) / A1(irow,k-1);
        temp2 = A2(irow,icol) / A2(irow,k-1);
        if (temp1 < minratios1(icol)){
          minratios1(icol) = temp1;
        }
        if (temp2 < minratios2(icol)){
          minratios2(icol) = temp2;
        }
      }
    }
  }
  return List::create(Named("pts1") = A1, Named("minratios1") = minratios1,
                      Named("pts2") = A2, Named("minratios2") = minratios2);
}

// Sample just one draw uniformly in the sub-simplex Delta_k(theta)
// [[Rcpp::export]]
NumericVector runif_piktheta_one_cpp(int k, const NumericVector & theta){
  RNGScope rngScope;
  int K = theta.length();
  GetRNGstate();
  NumericVector w = rexp(K, 1.0);
  PutRNGstate();
  w = w / sum(w);
  // at this point w is Dirichlet (1,..,1)
  for (int icol = 0; icol < K; icol ++){
    if (icol+1 != k){
      w(icol) = w(k-1) * theta(icol) + w(icol);
    }
  }
  w(k-1) *= theta(k-1);
  return w;
}


// Compute pdf of uniform distribution on Delta_k(theta) evaluated at x
// i.e. indicator(x in Delta_k(theta)) / volume(Delta_k(theta))
// [[Rcpp::export]]
double dunif_piktheta_cpp(const NumericVector & x_, int k, const NumericVector & theta){
  int K = theta.length();
  double pdf_ = 1./ theta(k-1); //  one over volume of subsimplex Delta_k(theta)
  for (int icol = 0; icol < K; icol ++){
    if (icol != k-1){
      if ((x_(icol) / x_(k-1)) < (theta(icol) / theta(k-1))){
        pdf_ = 0;
      }
    }
  }
  return pdf_;
}


// Maximum coupling of uniform distribution on sub-simplices, Delta_k(theta1) and Delta_k(theta2)
// [[Rcpp::export]]
List maxcoupling_runif_piktheta_cpp(int n, int k, const NumericVector & theta1, const NumericVector & theta2){
  RNGScope rngScope;
  int K = theta1.length();
  NumericMatrix A1(n, K);
  NumericMatrix A2(n, K);
  NumericVector x(K);
  NumericVector y(K);
  NumericVector u;
  double pdf1_x, pdf2_x, pdf1_y, pdf2_y;
  int ncoupled = 0;
  bool reject;
  double temp1, temp2;
  NumericVector minratios1(K);
  std::fill(minratios1.begin(), minratios1.end(), 0);
  NumericVector minratios2(K);
  std::fill(minratios2.begin(), minratios2.end(), 0);
  for (int irow = 0; irow < n; irow ++){
    // for each point ...
    x = runif_piktheta_one_cpp(k, theta1);
    pdf1_x = 1. / theta1(k-1);
    pdf2_x = dunif_piktheta_cpp(x, k, theta2);
    A1(irow,_) = clone(x);
    GetRNGstate();
    u = runif(1);
    PutRNGstate();
    if (u(0) < (pdf2_x / pdf1_x)){
      ncoupled ++;
      A2(irow,_) = clone(x);
    } else {
      reject = true;
      std::fill(y.begin(), y.end(), 0.);
      int iter = 0;
      y = runif_piktheta_one_cpp(k, theta2);
      while (iter < 100){
        iter ++;
        GetRNGstate();
        u = runif(1);
        PutRNGstate();      
        pdf2_y = 1. / theta2(k-1);
        pdf1_y = dunif_piktheta_cpp(y, k, theta1); 
        if (u(0) > (pdf1_y / pdf2_y)){
          iter = 100;
        } else {
          y = runif_piktheta_one_cpp(k, theta2);
        }
      }
      A2(irow,_) = clone(y);
    }
    // update minimum of A_n,ell/ A_n,k
    if (irow == 0){
      for (int icol = 0; icol < K; icol ++){
        minratios1(icol) = A1(0,icol) / A1(0,k-1);
        minratios2(icol) = A2(0,icol) / A2(0,k-1);
      }
    } else {
      for (int icol = 0; icol < K; icol ++){
        temp1 = A1(irow,icol) / A1(irow,k-1);
        temp2 = A2(irow,icol) / A2(irow,k-1);
        if (temp1 < minratios1(icol)){
          minratios1(icol) = temp1;
        }
        if (temp2 < minratios2(icol)){
          minratios2(icol) = temp2;
        }
      }
    }
  }
  return List::create(Named("pts1") = A1, Named("minratios1") = minratios1,
                      Named("pts2") = A2, Named("minratios2") = minratios2,
                      Named("ncoupled") = ncoupled);
}
// 
