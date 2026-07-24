# Normal prior on fixed selectivity parameters

Shared across the total fishery, retained fishery, and survey fixed
selectivity parameter priors in `SPoRC_rtmb.R` (the three "Selectivity
(Prior)" blocks) since all three prior tables and their corresponding
parameter arrays share the same `[region, par, block, sex, fleet]`
layout.

## Usage

``` r
get_selex_fixed_prior(selex_prior, fixed_sel_pars)
```

## Arguments

- selex_prior:

  Data frame with columns `region`, `par`, `block`, `sex`, `fleet`, `mu`
  (prior mean, natural scale), `sd` (prior SD, log scale) — one row per
  penalized parameter.

- fixed_sel_pars:

  Array `[region, par, block, sex, fleet]` of fixed selectivity
  parameters on the log scale.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `selex_prior`.
