# Validate a model-bin to observed-bin map

`AgeingError` and `LenBinMap` are the same operation on different axes:
an `[n_model_bins x n_obs_bins]` matrix the expected composition is
multiplied through so it lands on the bins the observations were
recorded on. The likelihood does not distinguish them, and neither does
this check.

## Usage

``` r
check_bin_map(x, n_model_bins, what, strict = TRUE, tol = 1e-08)
```

## Arguments

- x:

  The matrix to check.

- n_model_bins:

  Integer. Number of model bins, the required row count.

- what:

  Character. Argument name, used in messages.

- strict:

  Logical. `TRUE` (default) makes a bad row sum an error, `FALSE`
  reports it through
  [`collect_message`](https://chengmatt.github.io/SPoRC/dev/reference/collect_message.md).

- tol:

  Numeric. How far a row sum may sit from one before it is reported.

## Value

`x` invisibly, as a matrix.

## Details

A row is one model bin's share across the observed bins, so it sums to
one. A row of zeros is allowed and drops that model bin from the
observations, which is how observed bins starting above the first model
bin are expressed, as a shifted identity such as `diag(1, 10)[, 2:10]`.
A column of zeros is an observed bin nothing maps into, whose expected
proportion is a structural zero the composition likelihood cannot fit,
and a negative entry is fatal either way since nothing downstream can
read one.

The row-sum tolerance is the caller's choice. Published ageing error
matrices are rounded at source and come in with rows summing to 0.997 or
1.002, and the likelihood renormalizes after the multiply, so
`AgeingError` passes `tol = 0.05`. A length bin map is written by hand,
so `LenBinMap` keeps its `1e-8`. `strict` decides whether a row outside
the tolerance is fatal: `LenBinMap` rejects such a matrix outright,
while `AgeingError` reports it through the setup messages.
