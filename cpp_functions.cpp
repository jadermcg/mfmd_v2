#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector likelihood_t1C(List y, int w, NumericVector theta, NumericVector priori) {
  int L = y.size();
  NumericVector lhood(L);
  
  for(int i = 0; i < L; ++i) {
    NumericVector s = y[i];
    lhood[i] = priori[0];
    for(int j = 0; j < w; ++j) {
      lhood[i] *= theta[s[j] - 1];
    }
  }
  
  
  return lhood + 1e-50;
}

// [[Rcpp::export]]
NumericVector likelihood_t1C_v2(List Y, int w, NumericVector theta, NumericVector priori) {
  int N = Y.size();
  int L = static_cast<List>(Y[0]).size();
  
  NumericVector lhood(N*L);
  int z = 0;
  for(int k = 0; k < N; ++k) {
    List y = Y[k];
    for(int i = 0; i < L; ++i, ++z) {
      NumericVector seq = y[i];
      lhood[z] = priori[0];
      for(int j = 0; j < w; ++j) {
        lhood[z] *= theta[seq[j] - 1];
      }
    }
  }
  
  
  return lhood + 1e-50;
}

// [[Rcpp::export]]
NumericVector likelihood_t2C(List y, int w, NumericMatrix theta, NumericVector priori) {
  int L = y.size();
  NumericVector lhood(L);
  
  for(int i = 0; i < L; ++i) {
    NumericVector s = y[i];
    lhood[i] = priori[1];
    for(int j = 0; j < w; ++j) {
      lhood[i] *= theta(j, s[j] - 1);
    }
  }
  
  
  return lhood + 1e-50;
}

// [[Rcpp::export]]
NumericVector likelihood_t2C_v2(List Y, int w, NumericMatrix theta, NumericVector priori) {
  int N = Y.size();
  int L = static_cast<List>(Y[0]).size();
  
  NumericVector lhood(L*N);
  
  int z = 0;
  for(int k = 0; k < N; ++k) {
    List y = Y[k];
    for(int i = 0; i < L; ++i, ++z) {
      NumericVector seq = y[i];
      lhood[z] = priori[1];
      for(int j = 0; j < w; ++j) {
        lhood[z] *= theta(j, seq[j] - 1);
      }
    }
  }
  
  
  return lhood + 1e-50;
}

// [[Rcpp::export]]
NumericVector count_nucleotidesC(List Y, int w) {
  NumericVector counter{0,0,0,0};
  int n = Y.size();
  
  for(int i = 0; i < n; ++i) {
    NumericVector l  = Y[i];
    for(int j = 0; j < w; ++j) {
      ++counter[l[j] - 1];
    }
  }

  return counter;
}

// [[Rcpp::export]]
NumericMatrix w_mersC(NumericVector x, int w) {
  
  int m = x.size();
  int n = m - w + 1;
  NumericMatrix out{n,w};
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0, k = i; j < w; ++j, ++k) {
      out(i,j) = x[k];
    }
  }
  
  return out;
  
}

// [[Rcpp::export]]
NumericMatrix pfm_matrixC(NumericMatrix Y) {
  int n = Y.nrow();
  int m = Y.ncol();
  NumericMatrix pfm{m, 4};
  
  for(int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      if (Y(i,j) == 1) {
        ++pfm(j,0);
      } else if(Y(i,j) == 2) {
        ++pfm(j,1);
      } else if (Y(i,j) == 3) {
        ++pfm(j,2);
      } else if (Y(i,j) == 4) {
        ++pfm(j,3);
      }
    }
  }

  return pfm + 1;
  
}

// [[Rcpp::export]]
NumericMatrix update_pfmC(NumericMatrix pfm, NumericVector y) {
  int n = y.size();
  NumericMatrix new_pfm = clone(pfm);
  for(int i = 0; i < n; ++i) {
   ++new_pfm(i, y[i] - 1);
  }
  
  return new_pfm;
}

// [[Rcpp::export]]
NumericMatrix ppm_matrixC(NumericMatrix pfm) {
  int n = pfm.nrow();
  int m = pfm.ncol();
  double total = sum(pfm(0,_));
    
  return pfm / total;
  
}

// [[Rcpp::export]]
NumericMatrix pssm_matrixC(NumericVector theta_1, NumericMatrix theta_2) {
  int n = theta_2.nrow();
  int m = theta_2.ncol();
  NumericMatrix pssm{n, m};
  
  for(int i = 0; i < n; ++i) {
    pssm(i, _) = theta_2(i,_) * log(theta_2(i, _) / theta_1);
  }
  
  return pssm;
  
}

// [[Rcpp::export]]
NumericMatrix maximization_t2C(List Y, List rn2, int w) {
  NumericMatrix estimates(w, 4);
  int L = Y.size();
  int n = static_cast<List>(Y[0]).size();
  
  for(int k = 0; k < L; ++k) {
    List y = Y[k];
    NumericVector rn = rn2[k];
    for(int i = 0; i < n; ++i) {
      NumericVector seq = y[i];
      for(int j = 0; j < w; ++j) {
        switch( int(seq[j]) ) {
        case 1:
          estimates(j,0) += rn[i];
          break;
        case 2:
          estimates(j,1) += rn[i];
          break;
        case 3:
          estimates(j,2) += rn[i];
          break;
        case 4:
          estimates(j,3) += rn[i];
          break;
        }
      }
    }
  }
  
  
  
  return estimates;
}

// [[Rcpp::export]]
NumericMatrix maximization_t2C_v2(List Y, NumericVector rn2, int w) {
  NumericMatrix estimates(w, 4);
  int N = Y.size();
  int L = static_cast<List>(Y[0]).size();
  
  int z = 0;
  for(int k = 0; k < N; ++k) {
    List y = Y[k];
    for(int i = 0; i < L; ++i, ++z) {
      NumericVector seq = y[i];
      for(int j = 0; j < w; ++j) {
        switch( int(seq[j]) ) {
        case 1:
          estimates(j,0) += rn2[z];
          break;
        case 2:
          estimates(j,1) += rn2[z];
          break;
        case 3:
          estimates(j,2) += rn2[z];
          break;
        case 4:
          estimates(j,3) += rn2[z];
          break;
        }
      }
    }
  }
  return estimates;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
