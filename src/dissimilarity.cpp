// File: src/dissimilarity.cpp
// Description: Corrected C++ implementations of dissimilarity metrics for biodivMapR
//              Assumes A and B are abundance matrices. Binarizes for metrics expecting binary data.
// Author: Jean-Baptiste Feret
// Date: 2026-08-12

#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// Helper function to avoid division by zero
inline double safe_divide(double numerator, double denominator) {
  return (denominator == 0.0) ? 0.0 : numerator / denominator;
}

// Helper function to binarize a value (1 if > 0, else 0)
inline int binarize(double x) {
  return (x > 0) ? 1 : 0;
}

// --- Jaccard Dissimilarity (Binary) ---
// Formula: 1 - (intersection / union)
// [[Rcpp::export]]
NumericMatrix compute_jaccard_dissimilarity(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  int m = B.nrow();
  int p = A.ncol();

  NumericMatrix diss_matrix(n, m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      double intersection = 0.0;
      double union_set = 0.0;

      for (int k = 0; k < p; k++) {
        bool a_present = binarize(A(i, k));
        bool b_present = binarize(B(j, k));

        if (a_present || b_present) {
          union_set += 1.0;
        }
        if (a_present && b_present) {
          intersection += 1.0;
        }
      }
      diss_matrix(i, j) = 1.0 - safe_divide(intersection, union_set);
    }
  }
  return diss_matrix;
}

// --- Jaccard Turnover (Binary) ---
// Formula: (a + b) / (a + b + c)
// a: Species in A only, b: Species in B only, c: Species in both
// [[Rcpp::export]]
NumericMatrix compute_jaccard_turnover(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  int m = B.nrow();
  int p = A.ncol();

  NumericMatrix diss_matrix(n, m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      double a = 0.0; // Species in A only
      double b = 0.0; // Species in B only
      double c = 0.0; // Species in both

      for (int k = 0; k < p; k++) {
        bool a_present = binarize(A(i, k));
        bool b_present = binarize(B(j, k));

        if (a_present && !b_present) a += 1.0;
        if (!a_present && b_present) b += 1.0;
        if (a_present && b_present) c += 1.0;
      }
      diss_matrix(i, j) = safe_divide(a + b, a + b + c);
    }
  }
  return diss_matrix;
}

// --- Sorensen Dissimilarity (Binarized) ---
// Formula: 1 - (2 * intersection) / (sum(A_bin) + sum(B_bin))
// [[Rcpp::export]]
NumericMatrix compute_sorensen_dissimilarity(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  int m = B.nrow();
  int p = A.ncol();

  NumericMatrix diss_matrix(n, m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      double intersection = 0.0;
      double sum_A = 0.0;
      double sum_B = 0.0;

      for (int k = 0; k < p; k++) {
        // Binarize the data
        int a_bin = binarize(A(i, k));
        int b_bin = binarize(B(j, k));

        sum_A += a_bin;
        sum_B += b_bin;
        intersection += std::min(a_bin, b_bin);
      }
      diss_matrix(i, j) = 1.0 - safe_divide(2.0 * intersection, sum_A + sum_B);
    }
  }
  return diss_matrix;
}

// --- Simpson Dissimilarity (Quantitative) ---
// Formula: 1 - (intersection / min(sum(A), sum(B)))
// [[Rcpp::export]]
NumericMatrix compute_simpson_dissimilarity(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  int m = B.nrow();
  int p = A.ncol();

  NumericMatrix diss_matrix(n, m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      double intersection = 0.0;
      double sum_A = 0.0;
      double sum_B = 0.0;

      for (int k = 0; k < p; k++) {
        sum_A += A(i, k);
        sum_B += B(j, k);
        intersection += std::min(A(i, k), B(j, k));
      }
      diss_matrix(i, j) = 1.0 - safe_divide(intersection, std::min(sum_A, sum_B));
    }
  }
  return diss_matrix;
}

// --- Bray-Curtis Dissimilarity (Quantitative) ---
// Formula: sum(|A_ik - B_jk|) / sum(A_ik + B_jk)
// [[Rcpp::export]]
NumericMatrix compute_bray_curtis_dissimilarity(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  int m = B.nrow();
  int p = A.ncol();

  NumericMatrix diss_matrix(n, m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      double sum_abs_diff = 0.0;
      double sum_total = 0.0;

      for (int k = 0; k < p; k++) {
        sum_abs_diff += std::abs(A(i, k) - B(j, k));
        sum_total += A(i, k) + B(j, k);
      }
      diss_matrix(i, j) = safe_divide(sum_abs_diff, sum_total);
    }
  }
  return diss_matrix;
}

// --- Bray-Curtis Turnover (Quantitative) ---
// Formula: sum(|A_ik - B_jk|) / sum(max(A_ik, B_jk))
// [[Rcpp::export]]
NumericMatrix compute_bray_curtis_turnover(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  int m = B.nrow();
  int p = A.ncol();

  NumericMatrix diss_matrix(n, m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      double sum_abs_diff = 0.0;
      double sum_max = 0.0;

      for (int k = 0; k < p; k++) {
        sum_abs_diff += std::abs(A(i, k) - B(j, k));
        sum_max += std::max(A(i, k), B(j, k));
      }
      diss_matrix(i, j) = safe_divide(sum_abs_diff, sum_max);
    }
  }
  return diss_matrix;
}
