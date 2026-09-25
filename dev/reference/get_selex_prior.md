# Prior on selectivity, on the parameters or on realized values

Shared by the total fishery, retained fishery and survey selectivity
prior blocks in `SPoRC_rtmb.R`, since all three tables and their
parameter arrays are laid out over `[region, par, block, sex, fleet]`.
Each row is one prior, and its optional `type` column says what the row
constrains.

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
  `mu`, `sd` and optionally `type`, one row per prior. A `"par"` row
  (the default when the column is absent) is a lognormal prior on one
  fixed selectivity parameter,
  `dnorm(pars[region,par,block,sex,fleet], log(mu), sd)`, with `mu` on
  the natural scale and `sd` on the log scale. A `"value"` row is a
  normal prior on the realized selectivity at one bin,
  `dnorm(sel[bin], mu, sd)`, with both on the natural scale; `par` then
  names the bin on whichever grid the data source is parameterized over,
  and the value is read at the first model year of `block`, selectivity
  being constant within a block. A `"value"` row constrains a derived
  quantity rather than the parameters, which is the ADMB convention of
  pinning survey selectivity at a reference age near one, and expresses
  statements no set of independent parameter priors can, such as the
  rank-one ridge in (a50, slope) space a logistic curve's value at one
  age implies.

- fixed_sel_pars:

  Array `[region, par, block, sex, fleet]` of fixed selectivity
  parameters on the log scale, read by `"par"` rows.

- sel:

  Array `[pop, region, year, seas, age, sex, fleet]` of realized
  age-based selectivity, read by `"value"` rows at pop 1 and season 1,
  matching the smoothness penalties.

- sel_l:

  Array `[region, year, len, sex, fleet]` of realized length-based
  selectivity, read by `"value"` rows in place of `sel` when the data
  source is length-based.

- selex_type:

  Integer. `0` reads `sel`, `1` reads `sel_l`.

- sel_blocks:

  Integer array `[region, year, fleet]` mapping model years to
  selectivity blocks, resolving a `"value"` row's `block` to its first
  year.

## Value

Numeric scalar negative log-likelihood, summed over the rows of
`selex_prior`.
