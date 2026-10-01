# Compute Movement Process Error Log-Likelihood (Positive Scale)

The log likelihood of `move_devs`, a field over origin-destination pair,
population, year, season, age and sex (see
[`do_move_re_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_move_re_mapping.md)).
Each process error block of pairs is one separable Gaussian density,
`dseparable` with a factor per dim built from that dim's first level in
the block: `dautoreg` where the dim is `"ar1"`, `dmvnorm` with the
unstructured correlation where the dim is `"us"`, and a standard normal
otherwise. The sd in `PE_pars` is the conditional one, as `ln_sigmaNAA`
is for the numbers-at-age state, so for each `"ar1"` dim the field is
scaled by that sd divided by `sqrt(1 - rho^2)`, which gives the marginal
sd. With no correlation on any dim, each cell is penalized on its own,
and a cell left `NA` in `map_move_devs` is skipped – this is how a dsem
takes over a single series. With a correlation, a block is evaluated as
a whole: skipped when every cell is `NA`, refused when only some are.
Returns zero when the dsem holds the density (`move_dsem = 1`) or no
deviation is estimated.

## Usage

``` r
Get_move_PE_loglik(
  move_year_re,
  move_age_re,
  move_pop_re,
  move_seas_re,
  move_sex_re,
  PE_pars,
  move_pop_corr_pars,
  move_seas_corr_pars,
  move_sex_corr_pars,
  move_devs,
  map_move_devs,
  move_pairs,
  move_pe_block,
  move_pop_block,
  move_year_block,
  move_seas_block,
  move_age_block,
  move_sex_block,
  move_dsem
)
```

## Arguments

- move_year_re, move_age_re:

  Integer codes, `0` none, `1` iid, `2` ar1.

- move_pop_re, move_seas_re, move_sex_re:

  Integer codes, `0` none, `1` iid or blocks, `2` unstructured.

- PE_pars:

  Array `[n_regions x n_regions_to x 3]` of log conditional sd,
  unconstrained age correlation and unconstrained year correlation,
  shared within each process error block.

- move_pop_corr_pars, move_seas_corr_pars, move_sex_corr_pars:

  Unconstrained parameters of the unstructured correlations, read under
  `"us"`.

- move_devs:

  Array
  `[n_pop x n_regions x n_regions_to x n_years x n_seas x n_ages x n_sexes]`
  of movement deviations.

- map_move_devs:

  Integer mirror of the deviation map, `NA` where a cell is fixed or a
  dsem holds it.

- move_pairs:

  Integer matrix, one row per pair, giving its origin and destination.

- move_pe_block:

  Process error block of every pair.

- move_pop_block, move_year_block, move_seas_block, move_age_block,
  move_sex_block:

  Block of every level of the dim, `NA` where inactive.

- move_dsem:

  Integer flag; `1` when the dsem holds the density.

## Value

Scalar log likelihood (positive scale).
