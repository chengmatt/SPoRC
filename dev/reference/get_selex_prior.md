# Prior on selectivity, on the parameters or on realized values

Shared across the total fishery, retained fishery, and survey
"Selectivity (Prior)" blocks in `SPoRC_rtmb.R` since all three prior
tables and their corresponding parameter arrays share the same
`[region, par, block, sex, fleet]` layout. Each row of the table is one
prior, and its optional `type` column selects what the row constrains:

- `"par"` (the default when the column is absent):

  A lognormal prior on one fixed selectivity parameter,
  `dnorm(pars[region,par,block,sex,fleet], log(mu), sd)`, with `mu` on
  the natural scale and `sd` on the log scale.

- `"value"`:

  A normal prior on the realized selectivity value at one bin,
  `dnorm(sel[bin], mu, sd)`, with both hyperparameters on the natural
  scale. `par` instead names the bin, on the grid the stream's
  selectivity is parameterized on (ages or lengths per its selectivity
  type), and the value is read at the first model year of `block`
  (blocked and time-invariant selectivity are constant within a block).
  This is a constraint on a derived quantity rather than on the
  parameters (the ADMB rockfish convention of pinning survey selectivity
  at a reference age near one is its motivating case), so it can express
  statements no set of independent parameter priors can, e.g. the
  rank-one ridge in (a50, slope) space implied by constraining a
  logistic curve's value at one age.

## Usage

``` r
get_selex_prior(
  selex_prior,
  fixed_sel_pars,
  sel,
  sel_l,
  selex_type,
  sel_blocks
)
```

## Arguments

- selex_prior:

  Data frame with columns `region`, `par`, `block`, `sex`, `fleet`,
  `mu`, `sd`, and optionally `type` (`"par"`/`"value"`), one row per
  prior.

- fixed_sel_pars:

  Array `[region, par, block, sex, fleet]` of fixed selectivity
  parameters on the log scale, read by `"par"` rows.

- sel:

  Array `[pop, region, year, seas, age, sex, fleet]` of realized
  age-based selectivity, read by `"value"` rows at pop 1 and season 1,
  matching the smoothness penalties.

- sel_l:

  Array `[region, year, len, sex, fleet]` of realized length-based
  selectivity, read by `"value"` rows instead of `sel` when the stream
  is length-based.

- selex_type:

  Integer. `0` reads `sel`, `1` reads `sel_l`.

- sel_blocks:

  Integer array `[region, year, fleet]` mapping model years to
  selectivity blocks, used to resolve a `"value"` row's `block` to the
  first year in it.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `selex_prior`.
