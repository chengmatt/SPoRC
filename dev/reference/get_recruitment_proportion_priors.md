# Dirichlet/beta priors on recruitment apportionment

Combines the recruitment regional apportionment prior, seasonal
apportionment prior, and stray rate prior, since all three feed the
single `rec_prop_nLL` accumulator in `SPoRC_rtmb.R`. Called once from
the "Recruitment Proportions (Prior)" / "Stray Rates (Prior)" sections.

## Usage

``` r
get_recruitment_proportion_priors(
  use_rec_region_prop_prior,
  rec_region_prop_prior,
  rec_region_prop,
  use_rec_seas_prop_prior,
  use_fixed_rec_seas_prop,
  rec_seas_prop_prior,
  rec_seas_prop,
  rec_lag,
  spawn_seas,
  n_seas,
  use_stray_rate_prior,
  stray_rate_prior,
  stray_rate_pars
)
```

## Arguments

- use_rec_region_prop_prior:

  Integer (0/1) switch for the regional apportionment prior.

- rec_region_prop_prior:

  Data frame with columns `pop` and `alpha` (list column of Dirichlet
  concentration vectors).

- rec_region_prop:

  Array `[pop, region]` of recruitment regional apportionment.

- use_rec_seas_prop_prior, use_fixed_rec_seas_prop:

  Integer (0/1) switches; the seasonal apportionment prior is skipped
  when seasonal apportionment is fixed rather than estimated.

- rec_seas_prop_prior:

  Data frame with columns `pop` and `alpha` (list column of Dirichlet
  concentration vectors).

- rec_seas_prop:

  Array `[pop, season]` of recruitment seasonal apportionment.

- rec_lag, spawn_seas, n_seas:

  Integers controlling which seasons are structurally zero (age-0
  recruits before the spawning event) and so excluded from the seasonal
  Dirichlet prior.

- use_stray_rate_prior:

  Integer (0/1) switch for the stray rate prior.

- stray_rate_prior:

  Data frame with columns `pop`, `block`, `mu` (prior mean, natural
  scale), `sd` (prior SD, natural scale).

- stray_rate_pars:

  Array `[pop, block]` of stray rate parameters on the logit scale.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
three prior sources.
