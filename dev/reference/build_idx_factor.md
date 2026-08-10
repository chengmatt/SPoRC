# Precompute common-factor draw parameters for multivariate normal index fleets

Validates each multivariate normal fleet's covariance against its use
flags with
[`parse_idx_cov`](https://chengmatt.github.io/SPoRC/dev/reference/parse_idx_cov.md),
factor-decomposes it with
[`cov_to_factor`](https://chengmatt.github.io/SPoRC/dev/reference/cov_to_factor.md),
and records where each simulated cell sits in the covariance. Row `i` of
the covariance is the `i`-th cell with a use flag of 1 when scanning in
array order (region varies fastest, then year, then season), which is
the same order the estimation model collects the observation vector in,
so the simulated and fitted series line up.

## Usage

``` r
build_idx_factor(cov_list, like_type_vals, use_arr, n_fleets, what)
```

## Arguments

- cov_list:

  List with one element per fleet, each a covariance matrix or `NULL`.

- like_type_vals:

  Integer vector of index likelihood codes; only fleets coded `2` are
  decomposed.

- use_arr:

  Array `[region, year, season, fleet]` of use flags. Its year dimension
  may be shorter than the simulation when projection years extend past
  the data.

- n_fleets:

  Integer. Number of fleets.

- what:

  Character. Name used in error messages.

## Value

List with one element per fleet: `NULL` for non-mvn fleets, and
otherwise `d` and `lambda` from
[`cov_to_factor`](https://chengmatt.github.io/SPoRC/dev/reference/cov_to_factor.md),
a `row` lookup array `[region, year, season]` holding each used cell's
covariance row (`NA` elsewhere), and the means `d_mean` and
`lambda_mean` used for cells outside the covariance.
