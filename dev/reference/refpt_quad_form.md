# Delta method quadratic form

Evaluates \\d \Sigma d'\\. A precision matrix is applied by solving, so
the dense inverse is never built.

## Usage

``` r
refpt_quad_form(cov_obj, d)
```

## Arguments

- cov_obj:

  Output of
  [`refpt_par_cov`](https://chengmatt.github.io/SPoRC/dev/reference/refpt_par_cov.md).

- d:

  Numeric matrix `[n_quantity, n_par]` of sensitivities.

## Value

Numeric matrix `[n_quantity, n_quantity]`.
