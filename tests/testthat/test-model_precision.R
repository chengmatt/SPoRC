library(SPoRC)
library(Matrix)

test_that("Get_3d_precision works", {

  # Dimensions
  test_that("returns a sparse matrix of the correct dimension", {
    n_ages <- 4; n_yrs <- 5
    Q <- SPoRC:::Get_3d_precision(n_ages, n_yrs, 0.3, 0.2, 0.1, log(1), 1)
    expect_s4_class(Q, "sparseMatrix")
    expect_equal(dim(Q), c(n_ages * n_yrs, n_ages * n_yrs))
  })

  test_that("Q is square", {
    Q <- SPoRC:::Get_3d_precision(3, 6, 0.3, 0.2, 0.1, log(1), 1)
    expect_equal(nrow(Q), ncol(Q))
  })

  # Symmetry
  test_that("Q is symmetric — Var_Type 1 (conditional)", {
    Q <- SPoRC:::Get_3d_precision(4, 5, 0.4, 0.3, 0.2, log(2), 1)
    expect_equal(as.matrix(Q), t(as.matrix(Q)), tolerance = 1e-10)
  })

  test_that("Q is symmetric — Var_Type 0 (marginal / stationary)", {
    Q <- SPoRC:::Get_3d_precision(4, 5, 0.4, 0.3, 0.2, log(2), 0)
    expect_equal(as.matrix(Q), t(as.matrix(Q)), tolerance = 1e-10)
  })


  # Positive definiteness  (required for dgmrf use)
  test_that("Q is positive definite — Var_Type 1", {
    Q  <- SPoRC:::Get_3d_precision(4, 5, 0.4, 0.3, 0.2, log(1), 1)
    ev <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  })

  test_that("Q is positive definite — Var_Type 0", {
    Q  <- SPoRC:::Get_3d_precision(4, 5, 0.4, 0.3, 0.2, log(1), 0)
    ev <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  })


  # Variance parameterisation
  test_that("Var_Type 0: marginal variances are constant and equal exp(ln_var)", {
    ln_var <- log(3)
    Q      <- SPoRC:::Get_3d_precision(3, 4, 0.3, 0.2, 0.1, ln_var, 0)
    marg_vars <- diag(solve(as.matrix(Q)))
    expect_equal(marg_vars, rep(exp(ln_var), length(marg_vars)), tolerance = 1e-6)
  })

  test_that("larger ln_var inflates marginal variance (Var_Type 0)", {
    mv1 <- mean(diag(solve(as.matrix(SPoRC:::Get_3d_precision(3, 4, 0.3, 0.2, 0.1, log(1), 0)))))
    mv2 <- mean(diag(solve(as.matrix(SPoRC:::Get_3d_precision(3, 4, 0.3, 0.2, 0.1, log(4), 0)))))
    expect_gt(mv2, mv1)
  })


  # Zero partial correlations → diagonal Q (independent nodes)
  test_that("all zero pcorrs give diagonal Q — Var_Type 1", {
    Q   <- SPoRC:::Get_3d_precision(3, 4, 0, 0, 0, log(2), 1)
    off <- as.matrix(Q)[row(as.matrix(Q)) != col(as.matrix(Q))]
    expect_equal(off, rep(0, length(off)), tolerance = 1e-12)
  })

  test_that("all zero pcorrs give diagonal Q — Var_Type 0", {
    Q   <- SPoRC:::Get_3d_precision(3, 4, 0, 0, 0, log(2), 0)
    off <- as.matrix(Q)[row(as.matrix(Q)) != col(as.matrix(Q))]
    expect_equal(off, rep(0, length(off)), tolerance = 1e-12)
  })


  # Sparsity structure grows with added correlation dimensions
  test_that("adding age pcorr increases nnz relative to zero-pcorr baseline", {
    Q0 <- SPoRC:::Get_3d_precision(4, 5, 0,   0,   0,   0, 1)
    Qa <- SPoRC:::Get_3d_precision(4, 5, 0.4, 0,   0,   0, 1)
    expect_gt(Matrix::nnzero(Qa), Matrix::nnzero(Q0))
  })

  test_that("adding cohort pcorr further increases nnz", {
    Qa  <- SPoRC:::Get_3d_precision(4, 5, 0.4, 0.3, 0,   0, 1)
    Qac <- SPoRC:::Get_3d_precision(4, 5, 0.4, 0.3, 0.2, 0, 1)
    expect_gt(Matrix::nnzero(Qac), Matrix::nnzero(Qa))
  })

  # Axis wiring: each partial correlation reaches its own neighbour
  test_that("pcorr_age links age-adjacent nodes and pcorr_year year-adjacent nodes", {
    n_ages <- 3; n_yrs <- 4
    index  <- expand.grid(seq_len(n_ages), seq_len(n_yrs)) # node n is (age, year)

    # with a single partial correlation switched on every node has one
    # neighbour, so the only off-diagonal entries of Q are the pairs that
    # correlation links. Report the step each pair takes along the two axes.
    edge_steps <- function(pcorr_age, pcorr_year) {
      Q  <- as.matrix(SPoRC:::Get_3d_precision(n_ages, n_yrs, pcorr_age, pcorr_year, 0, 0, 1))
      nz <- which(abs(Q) > 1e-12 & upper.tri(Q), arr.ind = TRUE)
      unname(unique(cbind(abs(index[nz[,1],1] - index[nz[,2],1]),
                          abs(index[nz[,1],2] - index[nz[,2],2]))))
    }

    expect_equal(edge_steps(0.4, 0), matrix(c(1L, 0L), nrow = 1)) # one age apart, same year
    expect_equal(edge_steps(0, 0.4), matrix(c(0L, 1L), nrow = 1)) # same age, one year apart
  })

  # Asymmetry when n_ages != n_yrs and pcorr_age != pcorr_year
  test_that("swapping pcorr_age and pcorr_year changes Q when grid is non-square", {
    Q1 <- SPoRC:::Get_3d_precision(3, 7, 0.5, 0.2, 0.1, 0, 1)
    Q2 <- SPoRC:::Get_3d_precision(3, 7, 0.2, 0.5, 0.1, 0, 1)
    expect_false(isTRUE(all.equal(as.matrix(Q1), as.matrix(Q2))))
  })
})
