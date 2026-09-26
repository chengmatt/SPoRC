# Negative log density of the dsem grid

The grid is a multivariate normal with mean \\\mu\\ and precision \\Q =
(I - B)^{\top} V^{-1} (I - B)\\: \$\$-\ell = -\tfrac{1}{2}\log\|Q\| +
\tfrac{1}{2}(x - \mu)^{\top} Q (x - \mu) + \tfrac{n}{2}\log 2\pi,\$\$
evaluated by [`RTMB::dgmrf`](https://rdrr.io/pkg/RTMB/man/MVgauss.html).

## Usage

``` r
get_dsem_nLL(
  dsem_beta,
  ln_dsem_sd,
  x_grid,
  mu_grid,
  dsem_model,
  dsem_cells,
  delta0 = NULL,
  grid = NULL
)
```

## Arguments

- dsem_beta, ln_dsem_sd, dsem_model, dsem_cells, x_grid:

  As in
  [`get_dsem_matrices`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_matrices.md).

- mu_grid:

  Matrix `[year, series]` of series means: a covariate's mean, zero for
  a deviation series, and the log deterministic prediction for a numbers
  at age series. The density reads `x_grid - mu_grid`, so for numbers at
  age that difference is the innovation, which is stored nowhere.

- delta0:

  Optional numeric vector, one first-year offset per series.

- grid:

  Optional output of
  [`get_dsem_grid`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_grid.md)
  for this call, which the objective already holds, so the projection is
  not repeated.

## Value

Scalar negative log density.

## Details

A first-year offset \\\delta_0\\ shifts year one of each series and is
propagated through the paths by \\(I - B)^{-1}\\, which is how an AR1
started off its mean decays back at \\\rho^{t-1}\\.

A series whose sd is fixed at zero makes \\V\\ singular, so that \\Q\\
does not exist and the model is reduced rank. Its cells are solved out
of the others instead and the density is on what is left, with the
precision
[`get_dsem_Q_oo`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_Q_oo.md)
builds.
