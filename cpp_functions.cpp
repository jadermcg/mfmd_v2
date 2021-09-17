#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector likelihood_t1C_v3(List Y, int w, NumericVector theta_1, double priori) {
  int N = Y.size();
  int L = static_cast<List>(Y[0]).size();
  
  NumericVector lhood(N*L);
  int z = 0;
  for(int k = 0; k < N; ++k) {
    List y = Y[k];
    for(int i = 0; i < L; ++i, ++z) {
      NumericVector seq = y[i];
      lhood[z] = 1;
      for(int j = 0; j < w; ++j) {
        lhood[z] *= theta_1[seq[j] - 1];
      }
    }
  }
  
  return priori * lhood + 1e-50;
}

// [[Rcpp::export]]
NumericVector likelihood_t2C_v3(List Y, int w, NumericMatrix theta_2, double priori) {
  int N = Y.size();
  int L = static_cast<List>(Y[0]).size();
  
  NumericVector lhood(L*N);
  
  int z = 0;
  for(int k = 0; k < N; ++k) {
    List y = Y[k];
    for(int i = 0; i < L; ++i, ++z) {
      NumericVector seq = y[i];
      lhood[z] = 1;
      for(int j = 0; j < w; ++j) {
        lhood[z] *= theta_2(j, seq[j] - 1);
      }
    }
  }
  
  
  return priori * lhood + 1e-50;
}

// [[Rcpp::export]]
NumericVector sample_scoresC(List Y, NumericMatrix pssm, int w) {
  int N = Y.size();
  int L = static_cast<List>(Y[0]).size();
  NumericVector scores(N*L);
  
  int z = 0;
  for(int k = 0; k < N; ++k) {
    List y = Y[k];
    for(int i = 0; i < L; ++i, ++z) {
      NumericVector seq = y[i];
      for(int j = 0; j < w; ++j) {
        scores[z] += pssm(j, seq[j] - 1);
      }
    }
  }
  
  return scores;
} 

// [[Rcpp::export]]
NumericVector count_nucleotidesC(List Y, int w) {
  NumericVector counter(4);
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
NumericMatrix maximization_v2(List Y, NumericVector rn, int w) {
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
          estimates(j,0) += rn[z];
          break;
        case 2:
          estimates(j,1) += rn[z];
          break;
        case 3:
          estimates(j,2) += rn[z];
          break;
        case 4:
          estimates(j,3) += rn[z];
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
