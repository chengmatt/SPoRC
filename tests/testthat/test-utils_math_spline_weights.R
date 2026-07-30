library(SPoRC)
library(testthat)

test_that("Get_Natural_Cubic_Spline_Weights works", {

  test_that("matches stats::splinefun(method = 'natural')", {
    x_nodes <- c(0, 0.2, 0.5, 0.8, 1)
    y_nodes <- c(0.1, 0.5, 0.9, 0.6, 0.3)
    x_out <- seq(0, 1, length.out = 37)

    W <- Get_Natural_Cubic_Spline_Weights(x_nodes, x_out)
    mine <- as.vector(W %*% y_nodes)
    ref <- splinefun(x_nodes, y_nodes, method = "natural")(x_out)

    expect_equal(dim(W), c(length(x_out), length(x_nodes)))
    expect_equal(mine, ref, tolerance = 1e-8)
  })

  test_that("single node produces a constant curve", {
    x_out <- seq(0, 1, length.out = 11)
    W <- Get_Natural_Cubic_Spline_Weights(0.5, x_out)
    expect_equal(dim(W), c(length(x_out), 1))
    expect_true(all(as.vector(W %*% 3.2) == 3.2))
  })

  test_that("two nodes reduce to linear interpolation", {
    x_out <- seq(0, 1, length.out = 11)
    W <- Get_Natural_Cubic_Spline_Weights(c(0, 1), x_out)
    lin <- as.vector(W %*% c(2, 5))
    expect_equal(lin, 2 + 3 * x_out, tolerance = 1e-10)
  })

  test_that("evaluating exactly at the node positions recovers the node values", {
    x_nodes <- c(0, 0.3, 0.6, 1)
    y_nodes <- c(-1, 2, 0.5, 3)
    W <- Get_Natural_Cubic_Spline_Weights(x_nodes, x_nodes)
    expect_equal(as.vector(W %*% y_nodes), y_nodes, tolerance = 1e-10)
  })

  test_that("query points outside the node range are clamped, not extrapolated", {
    x_nodes <- c(0, 0.5, 1)
    y_nodes <- c(1, 2, 1)
    W_in  <- Get_Natural_Cubic_Spline_Weights(x_nodes, 1)
    W_out <- Get_Natural_Cubic_Spline_Weights(x_nodes, 1.5)
    expect_equal(as.vector(W_in %*% y_nodes), as.vector(W_out %*% y_nodes), tolerance = 1e-10)
  })

  test_that("errors on non-strictly-increasing node positions", {
    expect_error(Get_Natural_Cubic_Spline_Weights(c(0, 0.5, 0.5, 1), 1:2))
  })

})
