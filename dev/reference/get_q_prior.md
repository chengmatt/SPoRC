# Normal prior on log catchability

Shared across the fishery and survey catchability prior blocks in
`SPoRC_rtmb.R`, since both prior tables and their corresponding
catchability arrays share the same `[region, block, fleet]` layout.

## Usage

``` r
get_q_prior(q_prior, ln_q)
```

## Arguments

- q_prior:

  Data frame with columns `region`, `block`, `fleet`, `mu` (prior mean,
  natural scale), `sd` (prior SD, log scale) — one row per penalized
  parameter.

- ln_q:

  Array `[region, block, fleet]` of log catchability.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `q_prior`.
