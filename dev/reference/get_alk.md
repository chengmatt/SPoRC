# Binned age-length key

\\P(l \mid a)\\ on the length bins given by their lower edges: the
normal (or lognormal) CDF is taken at every lower edge, differenced, and
the tails below the first edge and above the last are accumulated into
the end bins, so every column sums to one.

## Usage

``` r
get_alk(len_lower, mu, sd, dist = 0)
```

## Arguments

- len_lower:

  Numeric vector of lower bin edges.

- mu:

  Mean length at age, possibly AD.

- sd:

  Spread of length at age, possibly AD.

- dist:

  Integer, 0 normal, 1 lognormal (`sd` on the log scale, mean corrected
  so the arithmetic mean stays `mu`).

## Value

Matrix `[n_lens x n_ages]`.
