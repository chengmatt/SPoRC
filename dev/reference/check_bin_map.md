# Validate a model-bin to observed-bin map

`AgeingError` and `LenBinMap` are the same operation on different axes:
an `[n_model_bins x n_obs_bins]` matrix that the expected composition is
multiplied through so it lands on the bins the observations were
recorded on. The likelihood does not distinguish them, and neither does
this check, so a mistake in either one is reported the same way.

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
observations entirely, which is how observed bins that start above the
first model bin are expressed (a shifted identity such as
`diag(1, 10)[, 2:10]`).

The row-sum tolerance is a caller's choice. Published ageing error
matrices are rounded at source, and real ones come in with rows summing
to 0.997 or 1.002; the likelihood renormalizes the expectation after the
multiply, so a row off by that much reweights nothing, and `AgeingError`
passes `tol = 0.05`. A length bin map is written by hand rather than
read from a rounded table, so `LenBinMap` keeps the `1e-8` it has always
been held to. Only a row off by more than `tol` is reported, since that
means the matrix is not the map its author thought it was.

`strict` decides whether that is fatal. `LenBinMap` has always rejected
such a matrix outright and keeps doing so. `AgeingError` has not been
checked before, so a bad row is reported through the setup messages
rather than stopping a model that ran yesterday.

A column of zeros is an observed bin nothing maps into, whose expected
proportion is a structural zero the composition likelihood cannot fit.
It follows `strict` for the same reason the row sums do. A negative
entry is fatal either way, since nothing downstream can interpret one.
