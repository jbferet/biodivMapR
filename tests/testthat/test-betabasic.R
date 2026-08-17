test_that("Beta diversity for simple cases", {
    # Create a presence-absence matrix with identical sites
  A <- matrix(c(1, 1, 0, 0), nrow = 1, byrow = TRUE)
  B1 <- matrix(c(1, 1, 0, 0), nrow = 1, byrow = TRUE)
  colnames(A) <- c("Sp1", "Sp2", "Sp3", "Sp4")
  bet <- biodivMapR::compute_dissimilarity(A = A, B = B1,
                                           beta_metrics = c('bray', 'brayturn', 'simpson_diss',
                                                            'sorensen', 'jaccard', 'jaccardturn'))
  # Test for Jaccard, Sorensen, Bray-Curtis, etc.
  testthat::expect_equal(c(bet$bray), 0)
  testthat::expect_equal(c(bet$brayturn), 0)
  testthat::expect_equal(c(bet$simpson_diss), 0)
  testthat::expect_equal(c(bet$sorensen), 0)
  testthat::expect_equal(c(bet$jaccard), 0)
  testthat::expect_equal(c(bet$jaccardturn), 0)

  B2 <- matrix(c(0, 0, 1, 1), nrow = 1, byrow = TRUE)
  bet <- biodivMapR::compute_dissimilarity(A = A, B = B2,
                                           beta_metrics = c('bray', 'brayturn', 'simpson_diss',
                                                            'sorensen', 'jaccard', 'jaccardturn'))
  # Test for Jaccard, Sorensen, Bray-Curtis, etc.
  testthat::expect_equal(c(bet$bray), 1)
  testthat::expect_equal(c(bet$brayturn), 1)
  testthat::expect_equal(c(bet$simpson_diss), 1)
  testthat::expect_equal(c(bet$sorensen), 1)
  testthat::expect_equal(c(bet$jaccard), 1)
  testthat::expect_equal(c(bet$jaccardturn), 1)

  B3 <- matrix(c(0, 1, 1, 0), nrow = 1, byrow = TRUE)
  bet <- biodivMapR::compute_dissimilarity(A = A, B = B3,
                                           beta_metrics = c('bray', 'brayturn', 'simpson_diss',
                                                            'sorensen', 'jaccard', 'jaccardturn'))
  # Test for Jaccard, Sorensen, Bray-Curtis, etc.
  testthat::expect_equal(c(bet$simpson_diss), 0.5, tolerance = 0.01)
  testthat::expect_equal(c(bet$sorensen), 0.5, tolerance = 0.01)
  testthat::expect_equal(c(bet$jaccard), 2/3, tolerance = 0.01)
  testthat::expect_equal(c(bet$jaccardturn), 2/3, tolerance = 0.01)

  A2 <- matrix(c(10, 20, 0), nrow = 1, byrow = TRUE)
  B4 <- matrix(c(0, 15, 25), nrow = 1, byrow = TRUE)
  bet <- biodivMapR::compute_dissimilarity(A = A2, B = B4,
                                           beta_metrics = c('bray', 'brayturn', 'simpson_diss',
                                                            'sorensen', 'jaccard', 'jaccardturn'))
  # Test for Jaccard, Sorensen, Bray-Curtis, etc.
  testthat::expect_equal(c(bet$bray), 1-(2*15)/70, tolerance = 0.01)
  testthat::expect_equal(c(bet$brayturn), 40/55, tolerance = 0.01)
})
