# Conditional distribution of unknown cells in a Gaussian Markov random field

For \\x \sim N(\mu, Q^{-1})\\ split into unknown cells \\u\\ and known
cells \\k\\, \\x_u \mid x_k \sim N(\mu_u - Q\_{uu}^{-1} Q\_{uk} (x_k -
\mu_k), Q\_{uu}^{-1})\\. Only the unknown block of the precision is
factored.

## Usage

``` r
get_dsem_conditional(Q, mu_cell, known_cell, x_known)
```

## Arguments

- Q:

  Sparse precision matrix.

- mu_cell:

  Mean of every cell.

- known_cell:

  Indices of the known cells.

- x_known:

  Values of the known cells, in the order of `known_cell`.

## Value

List with `unknown_cell`, `cond_mean` and `chol_uu`, the sparse Cholesky
factor of \\Q\_{uu}\\.
