// File: src/compute_alpha_diversity_metrics.cpp
// Description: C++ implementations of alpha diversity metrics for biodivMapR
//              abundance matrix : species = columns, sites/samples = rows
//              compute_richness boolean
//              compute_shannon boolean
//              compute_simpson boolean
//              compute_hill boolean
//              q numeric
// Author: Jean-Baptiste Feret
// Date: 2026-08-12

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <numeric>

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_alpha_diversity_metrics(
    Rcpp::NumericMatrix abundances,
    bool compute_richness = true,
    bool compute_shannon = true,
    bool compute_simpson = true,
    bool compute_hill = true,
    double q = 1.0  // Default q value for Hill number (e.g., q=1 for Shannon-based Hill number)
) {
  int n_rows = abundances.nrow();

  // Count the number of columns needed in the result
  int n_cols = 0;
  if (compute_richness) n_cols++;
  if (compute_shannon) n_cols++;
  if (compute_simpson) n_cols++;
  if (compute_hill) n_cols++;
  // {
  //   // For Hill numbers: q=0,1,2 are always included if compute_hill=true
  //   // Custom q_values are added if provided
  //   n_cols += 3; // q=0,1,2
  //   for (double q : q_values) {
  //     if (q != 0 && q != 1 && q != 2) {
  //       n_cols++;
  //     }
  //   }
  // }

  Rcpp::NumericMatrix results(n_rows, n_cols);
  Rcpp::CharacterVector colnames;

  // Add column names based on selected metrics
  if (compute_richness) colnames.push_back("Richness");
  if (compute_shannon) colnames.push_back("Shannon");
  if (compute_simpson) colnames.push_back("Simpson");
  if (compute_hill) colnames.push_back("Hill");
  //   {
  //   colnames.push_back("Hill");
  //   // colnames.push_back("Hill_q1");
  //   // colnames.push_back("Hill_q2");
  //   for (double q : q_values) {
  //     if (q != 0 && q != 1 && q != 2) {
  //       colnames.push_back("Hill_q" + std::to_string(static_cast<int>(q)));
  //     }
  //   }
  // }

  // Set column names for the result matrix
  results.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);

  for (int i = 0; i < n_rows; i++) {
    Rcpp::NumericVector row = abundances(i, Rcpp::_);
    double total = std::accumulate(row.begin(), row.end(), 0.0);

    int col_index = 0;

    // Compute Richness
    if (compute_richness) {
      int richness = 0;
      for (double abundance : row) {
        if (abundance > 0) {
          richness++;
        }
      }
      results(i, col_index++) = static_cast<double>(richness);
    }

    // Compute Shannon Index
    if (compute_shannon) {
      if (total == 0) {
        results(i, col_index++) = 0.0;
      } else {
        double shannon = 0.0;
        for (double abundance : row) {
          if (abundance > 0) {
            double proportion = abundance / total;
            shannon -= proportion * std::log(proportion);
          }
        }
        results(i, col_index++) = shannon;
      }
    }

    // Compute Simpson Index
    if (compute_simpson) {
      if (total == 0) {
        results(i, col_index++) = 0.0;
      } else {
        double sum_squares = 0.0;
        for (double abundance : row) {
          if (abundance > 0) {
            double proportion = abundance / total;
            sum_squares += proportion * proportion;
          }
        }
        results(i, col_index++) = 1.0 - sum_squares;
      }
    }

    // Compute Hill Numbers
    if (compute_hill) {
      // Precompute richness, shannon, and sum_squares for Hill numbers
      int richness = 0;
      double shannon = 0.0;
      double sum_squares = 0.0;
      for (double abundance : row) {
        if (abundance > 0) {
          richness++;
          double proportion = abundance / total;
          shannon -= proportion * std::log(proportion);
          sum_squares += proportion * proportion;
        }
      }

      // Hill q=0 (richness)
      if (q == 0) {
        results(i, col_index++) = static_cast<double>(richness);
      }

      // Hill q=1 (exp(Shannon))
      if (q == 1) {
        results(i, col_index++) = (total > 0) ? std::exp(shannon) : 0.0;
      }

      // Hill q=2 (1 / sum_squares)
      if (q == 2) {
        results(i, col_index++) = (total > 0) ? 1.0 / sum_squares : 0.0;
      }

      // Custom Hill numbers for additional q values
      if (q == 0 || q == 1 || q == 2) continue; // Skip if already computed
      if (total == 0) {
        results(i, col_index++) = 0.0;
        continue;
      }
      double sum_pow = 0.0;
      for (double abundance : row) {
        if (abundance > 0) {
          double proportion = abundance / total;
          sum_pow += std::pow(proportion, q);
        }
      }
      results(i, col_index++) = std::pow(sum_pow, 1.0 / (1.0 - q));
    }
  }

  return results;}
