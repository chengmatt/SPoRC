# Simulate age or length compositions

Draws observed composition samples (by age or length) for a single
region–year–fleet–season–simulation cell, supporting multinomial,
Dirichlet-multinomial, and logistic-normal likelihoods. Ageing error is
optionally applied post-draw. Three composition aggregation structures
are handled: sex-split (`comp_type = 1`), joint across sexes
(`comp_type = 2`), and spatially aggregated across all regions
(`comp_type = 0`). The sentinel value `comp_type = 999` or
`comp_like = 999` causes the function to return `Obs` unchanged.

## Usage

``` r
simulate_comps(
  r,
  y,
  f,
  seas,
  sim,
  Exp,
  ISS = NULL,
  AgeingError,
  comp_like = NULL,
  ln_theta = NULL,
  corr_pars = NULL,
  ln_theta_agg = NULL,
  corr_pars_agg = NULL,
  comp_type = NULL,
  n_sexes,
  n_pop = NULL,
  n_regions,
  n_cat,
  Obs,
  pop_specific = FALSE,
  ISS_pop = NULL,
  pop_comp_like = NULL,
  pop_comp_type = NULL,
  ln_pop_theta = NULL,
  pop_corr_pars = NULL,
  ln_pop_theta_agg = NULL,
  pop_corr_pars_agg = NULL,
  age_or_len = 0
)
```

## Arguments

- r:

  Integer. Region index.

- y:

  Integer. Year index.

- f:

  Integer. Fleet index (fishery or survey).

- seas:

  Integer. Season index.

- sim:

  Integer. Simulation replicate index.

- Exp:

  Array. Expected compositions
  `[n_pop × n_regions × n_yrs × n_seas × n_cat × n_sexes × n_fleets × n_sims]`.

- ISS:

  Array. Integer sample sizes
  `[n_regions × n_yrs × n_seas × n_sexes × n_fleets × n_sims]`. Used
  when `pop_specific = FALSE`.

- AgeingError:

  Array. Ageing error transition matrices
  `[n_yrs × n_obs_ages × n_ages × n_sims]`. Ignored when
  `age_or_len = 1`.

- comp_like:

  Integer vector `[n_fleets]`. Likelihood type per fleet: `0` =
  multinomial, `1` = Dirichlet-multinomial, `2`–`4` = logistic-normal
  variants.

- ln_theta:

  Array. Log overdispersion or log-variance parameters
  `[n_regions × n_sexes × n_fleets]`. Used when `pop_specific = FALSE`.

- corr_pars:

  Array. Correlation parameters for logistic-normal likelihoods
  `[n_regions × n_sexes × n_fleets × n_corr_pars]`.

- ln_theta_agg:

  Numeric vector `[n_fleets]`. Log overdispersion for spatially
  aggregated compositions (`comp_type = 0`).

- corr_pars_agg:

  Numeric vector `[n_fleets]`. Correlation parameter(s) for aggregated
  logistic-normal compositions.

- comp_type:

  Integer matrix `[n_yrs × n_fleets]`. Aggregation structure: `0` =
  aggregated across regions, `1` = split by sex, `2` = joint across
  sexes, `999` = no data (skip).

- n_sexes:

  Integer. Number of sexes.

- n_pop:

  Integer. Number of populations.

- n_regions:

  Integer. Number of regions.

- n_cat:

  Integer. Number of composition categories (ages or lengths).

- Obs:

  Array. Observed compositions container with the same dimensions as
  `Exp`. Simulated values are written in-place.

- pop_specific:

  Logical. If `TRUE`, simulate compositions separately for each
  population using population-specific inputs.

- ISS_pop:

  Array. Population-specific sample sizes
  `[n_pop × n_regions × n_yrs × n_seas × n_sexes × n_fleets × n_sims]`.
  Used when `pop_specific = TRUE`.

- pop_comp_like:

  Integer vector `[n_fleets]`. Likelihood type per fleet for
  population-specific compositions.

- pop_comp_type:

  Integer matrix `[n_yrs × n_fleets]`. Aggregation structure for
  population-specific compositions.

- ln_pop_theta:

  Array. Log overdispersion parameters
  `[n_pop × n_regions × n_sexes × n_fleets]`.

- pop_corr_pars:

  Array. Correlation parameters for logistic-normal likelihoods
  `[n_pop × n_regions × n_sexes × n_fleets × n_corr_pars]`.

- ln_pop_theta_agg:

  Numeric array `[n_pop × n_fleets]`. Log overdispersion for
  population-specific aggregated compositions.

- pop_corr_pars_agg:

  Numeric array `[n_pop × n_fleets]`. Correlation parameter(s) for
  population-specific aggregated logistic-normal compositions.

- age_or_len:

  Integer. Indicator for composition type: `0` = age compositions (apply
  ageing error), `1` = length compositions (no ageing error).

## Value

The `Obs` array with simulated composition draws filled in at the
appropriate slice. When `pop_specific = FALSE`, values are written to
`[r, y, seas, , , f, sim]`; when `pop_specific = TRUE`, values are
written to `[p, r, y, seas, , , f, sim]`. All other slices are
unchanged.

## Details

When `pop_specific = TRUE`, compositions are simulated separately for
each population, extending all relevant inputs (e.g., sample size,
dispersion, and correlation parameters) to include a population
dimension. In this case, aggregation across regions (`comp_type = 0`) is
performed within population, and results are written to `Obs[p, ...]`.

For joint compositions (`comp_type = 2`), the Kronecker product
`diag(n_sexes) ⊗ AgeingError` is used to apply ageing error across the
combined age–sex vector. For aggregated compositions (`comp_type = 0`),
the draw is only executed when `r == n_regions` (i.e., on the final
region pass), and uses region- and sex-marginalised expected
proportions.
