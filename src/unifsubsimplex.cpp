// # dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
// # Copyright (C) 2019 Pierre E. Jacob
// # 
// # This program is free software: you can redistribute it and/or modify
// # it under the terms of the GNU General Public License as published by
// # the Free Software Foundation, either version 3 of the License, or
// # (at your option) any later version.
// # 
// # This program is distributed in the hope that it will be useful,
// # but WITHOUT ANY WARRANTY; without even the implied warranty of
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// # GNU General Public License for more details.
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <RcppEigen.h>
using namespace Rcpp;

// the following draws U_n uniformly in pi_k(theta), for all n
// and also computes the minimum of U_n,ell / U_n,k  for all n,ell

// [[Rcpp::export]]
List runif_piktheta_cpp(int n, int k, NumericVector theta){
  RNGScope rngScope;
  int K = theta.length();
  GetRNGstate();
  NumericVector w = rexp(n * K, 1.0);
  PutRNGstate();
  NumericMatrix U(n, K);
  double sumrow = 0.;
  double temp = 0.;
  NumericVector minratios(K);
  std::fill(minratios.begin(), minratios.end(), 0);
  for (int irow = 0; irow < n; irow ++){
    sumrow = 0.;
    for (int icol = 0; icol < K; icol ++){
      U(irow,icol) = w(irow*K + icol);
      sumrow += U(irow,icol);
    }
    for (int icol = 0; icol < K; icol ++){
      U(irow,icol) /= sumrow;
      // at this point U(irow,icol) is Dirichlet (1,..,1)
    }
    for (int icol = 0; icol < K; icol ++){
      if (icol+1 != k){
        U(irow,icol) = U(irow,k-1) * theta(icol) + U(irow,icol);
      }
    }
    U(irow,k-1) *= theta(k-1);
    // at this point we have U_n
    // update minimum of U_n,ell/ U_n,k
    if (irow == 0){
      for (int icol = 0; icol < K; icol ++){
        minratios(icol) = U(0,icol) / U(0,k-1);
      }
    } else {
      for (int icol = 0; icol < K; icol ++){
        temp = U(irow,icol) / U(irow,k-1);
        if (temp < minratios(icol)){
          minratios(icol) = temp;
        }
      }
    }
  }
  return List::create(Named("pts") = U, Named("minratios") = minratios);
}

// Coupled version, with common random numbers
// [[Rcpp::export]]
List crng_runif_piktheta_cpp(int n, int k, NumericVector & theta1, NumericVector & theta2){
  RNGScope rngScope;
  int K = theta1.length();
  GetRNGstate();
  NumericVector w = rexp(n * K, 1.0);
  PutRNGstate();
  NumericMatrix U1(n, K);
  NumericMatrix U2(n, K);
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
      U1(irow,icol) = w(irow*K + icol);
      sumrow1 += U1(irow,icol);
      U2(irow,icol) = w(irow*K + icol);
      sumrow2 += U2(irow,icol);
    }
    for (int icol = 0; icol < K; icol ++){
      U1(irow,icol) /= sumrow1;
      U2(irow,icol) /= sumrow2;
      // at this point U(irow,icol) is Dirichlet (1,..,1)
    }
    for (int icol = 0; icol < K; icol ++){
      if (icol+1 != k){
        U1(irow,icol) = U1(irow,k-1) * theta1(icol) + U1(irow,icol);
        U2(irow,icol) = U2(irow,k-1) * theta2(icol) + U2(irow,icol);
      }
    }
    U1(irow,k-1) *= theta1(k-1);
    U2(irow,k-1) *= theta2(k-1);
    // at this point we have U_n
    // update minimum of U_n,ell/ U_n,k
    if (irow == 0){
      for (int icol = 0; icol < K; icol ++){
        minratios1(icol) = U1(0,icol) / U1(0,k-1);
        minratios2(icol) = U2(0,icol) / U2(0,k-1);
      }
    } else {
      for (int icol = 0; icol < K; icol ++){
        temp1 = U1(irow,icol) / U1(irow,k-1);
        temp2 = U2(irow,icol) / U2(irow,k-1);
        if (temp1 < minratios1(icol)){
          minratios1(icol) = temp1;
        }
        if (temp2 < minratios2(icol)){
          minratios2(icol) = temp2;
        }
      }
    }
  }
  return List::create(Named("pts1") = U1, Named("minratios1") = minratios1,
                      Named("pts2") = U2, Named("minratios2") = minratios2);
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
// the 'multiplier' allows for a sub-maximal coupling
// but with a cost that has a bounded variance no matter what (following discussion with Anthony Lee)
// [[Rcpp::export]]
List maxcoupling_runif_piktheta_cpp(int n, int k, const NumericVector & theta1, const NumericVector & theta2){
  double multiplier = 0.95;
  RNGScope rngScope;
  int K = theta1.length();
  NumericMatrix U1(n, K);
  NumericMatrix U2(n, K);
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
  LogicalVector equal(n);
  for (int irow = 0; irow < n; irow ++){
    // for each point ...
    x = runif_piktheta_one_cpp(k, theta1);
    pdf1_x = 1. / theta1(k-1);
    pdf2_x = dunif_piktheta_cpp(x, k, theta2);
    U1(irow,_) = clone(x);
    GetRNGstate();
    u = runif(1);
    PutRNGstate();
    if (u(0) < (multiplier * pdf2_x / pdf1_x)){
      ncoupled ++;
      U2(irow,_) = clone(x);
      equal(irow) = true;
    } else {
      reject = true;
      std::fill(y.begin(), y.end(), 0.);
      int iter = 0;
      y = runif_piktheta_one_cpp(k, theta2);
      bool finish = false;
      while (!finish){
        iter ++;
        GetRNGstate();
        u = runif(1);
        PutRNGstate();      
        pdf2_y = 1. / theta2(k-1);
        pdf1_y = dunif_piktheta_cpp(y, k, theta1); 
        double acceptratio = (pdf1_y / pdf2_y);
        if (multiplier < acceptratio){
          acceptratio = multiplier;
        }
        if (u(0) > acceptratio){
          // finish
          finish = true;
        } else {
          // try again
          y = runif_piktheta_one_cpp(k, theta2);
        }
      }
      U2(irow,_) = clone(y);
      equal(irow) = false;
    }
    // update minimum of U_n,ell/ U_n,k
    if (irow == 0){
      for (int icol = 0; icol < K; icol ++){
        minratios1(icol) = U1(0,icol) / U1(0,k-1);
        minratios2(icol) = U2(0,icol) / U2(0,k-1);
      }
    } else {
      for (int icol = 0; icol < K; icol ++){
        temp1 = U1(irow,icol) / U1(irow,k-1);
        temp2 = U2(irow,icol) / U2(irow,k-1);
        if (temp1 < minratios1(icol)){
          minratios1(icol) = temp1;
        }
        if (temp2 < minratios2(icol)){
          minratios2(icol) = temp2;
        }
      }
    }
  }
  return List::create(Named("pts1") = U1, Named("minratios1") = minratios1,
                      Named("pts2") = U2, Named("minratios2") = minratios2,
                      Named("ncoupled") = ncoupled, Named("equal") = equal);
}
// 
