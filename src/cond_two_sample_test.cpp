#include <Rcpp.h>
using namespace Rcpp;

// Function to perform comparison
// [[Rcpp::export]]
NumericVector comp_cpp(NumericVector x) {
  int n = x.size() - 1;
  NumericVector y(n, x[0]);
  NumericVector x1 = x[Range(1, n)];
  for (int i = 0; i < n; ++i) {
    if (x1[i] > x[0]) {
      y[i] = x1[i];
    }
  }
  return y;
}

// Function to compute the Euclidean distance matrix
// [[Rcpp::export]]
NumericMatrix distMatrix_cpp(NumericMatrix V) {
  int m = V.nrow();
  NumericMatrix D(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      double sumSq = 0.0;
      for (int k = 0; k < V.ncol(); ++k) {
        double diff = V(i, k) - V(j, k);
        sumSq += diff * diff;
      }
      D(i, j) = sqrt(sumSq);
    }
  }
  return D;
}

// Function to combine matrices by row
// [[Rcpp::export]]
NumericMatrix combineRows_cpp(NumericMatrix A, NumericMatrix B) {
  int rowsA = A.nrow(), rowsB = B.nrow(), cols = A.ncol();
  NumericMatrix result(rowsA + rowsB, cols);
  for (int i = 0; i < rowsA; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(i, j) = A(i, j);
    }
  }
  for (int i = 0; i < rowsB; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(rowsA + i, j) = B(i, j);
    }
  }
  return result;
}

// Function to compute the statistic
// [[Rcpp::export]]
double Stat_cpp(NumericMatrix D) {
  int m = D.nrow();
  int n = m/2;
  // int n = X.nrow();
  // NumericMatrix X1 = cbind(X, Z);
  // NumericMatrix Y1 = cbind(Y, Z);
  // NumericMatrix V = combineRows_cpp(X1, Y1);
  // NumericMatrix D = distMatrix_cpp(V);
  double s = median(D);
  
  NumericVector T(n);
  // Loop over i = 0 to n-1 (corresponding to R's 1:n)
  for (int i = 0; i < n; ++i) {
    // Allocate matrices with n+1 columns (first column from a single vector,
    // then n columns from a block)  
    NumericMatrix M1(m, n + 1),
    M2(m, n + 1),
    M3(m, n + 1);
    for (int j = 0; j < m; ++j) {
      // For M1: first column is D(j, i), then columns 1..n from D(j, 0..n-1)
      M1(j, 0) = D(j, i) / s;
      for (int k = 1; k <= n; ++k) {
        M1(j, k) = D(j, k - 1) / s;
      }
      // For M2: first column is D(j, i), then columns 1..n come from D(j, n..2*n-1)
      M2(j, 0) = D(j, i) / s;
      for (int k = 1; k <= n; ++k) {
        M2(j, k) = D(j, (k - 1) + n) / s;
      }
      // For M3: first column is D(j, i+n), then columns 1..n again from D(j, n..2*n-1)
      M3(j, 0) = D(j, i + n) / s;
      for (int k = 1; k <= n; ++k) {
        M3(j, k) = D(j, (k - 1) + n) / s;
      }
    }
    // Now compute the row-wise sums:
    double sumM1 = 0.0, sumM2 = 0.0, sumM3 = 0.0;
    // Loop over each row of the M matrices.
    for (int j = 0; j < m; j++) {
      // Extract row j as a NumericVector:
      NumericVector row1 = M1(j, _);
      NumericVector row2 = M2(j, _);
      NumericVector row3 = M3(j, _);
      
      // Apply comp to each row:
      NumericVector comp1 = comp_cpp(row1);
      NumericVector comp2 = comp_cpp(row2);
      NumericVector comp3 = comp_cpp(row3);
      
      // Now, apply F (here, the normal CDF using R::pnorm) elementwise and sum.
      // Here we loop over the elements of the comp output.
      for (int l = 0; l < comp1.size(); l++) {
        sumM1 += R::pexp(comp1[l], 1.0, TRUE, FALSE);
        sumM2 += R::pexp(comp2[l], 1.0, TRUE, FALSE);
        sumM3 += R::pexp(comp3[l], 1.0, TRUE, FALSE);
      }
    }
    // For index i, assign:
    T[i] = 2 * sumM2 - sumM1 - sumM3;
  }
  return sum(T) / (n * n * n);
}

// Function to compute the test
// [[Rcpp::export]]
double Test_cpp(NumericMatrix X, NumericMatrix Y, NumericMatrix Z, int R) {
  int n = X.nrow();
  NumericMatrix X1 = cbind(X, Z);
  NumericMatrix Y1 = cbind(Y, Z);
  NumericMatrix V = combineRows_cpp(X1, Y1);
  NumericMatrix D = distMatrix_cpp(V);
  double S = Stat_cpp(D);
  NumericVector S1(R);
  
  for (int r = 0; r < R; ++r) {
    IntegerVector v = as<IntegerVector>(Rcpp::rbinom(n, 1, 0.5));
    std::vector<int> index;
    for (int i = 0; i < v.size(); ++i) {
      if (v[i] == 1) {
        index.push_back(i);
      }
    }
    IntegerVector pi(2 * n);
    for (int i = 0; i < n; ++i) {
      pi[i] = i;
      pi[n + i] = n + i;
    }
    for (int i = 0; i < index.size(); ++i) {
      std::swap(pi[index[i]], pi[n + index[i]]);
    }
    NumericMatrix D1(D.nrow(), D.ncol());
    for (int i = 0; i < D.nrow(); ++i) {
      for (int j = 0; j < D.ncol(); ++j) {
        D1(i, j) = D(pi[i], pi[j]);
      }
    }
    S1[r] = Stat_cpp(D1);
  }
  
  return (float) (1 + sum(S1 > S) + sum(S1 == S)) / (1 + R);
}
