# Simulate age or length compositions

Draws composition samples for one region, year, fleet, season and
replicate under the multinomial, Dirichlet-multinomial or
logistic-normal likelihoods, applying ageing error after the draw for
age compositions. A `comp_type` or `comp_like` of `999` returns `Obs`
unchanged.

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

- r, y, f, seas, sim:

  Region, year, fleet, season and replicate indices.

- Exp:

  Expected compositions
  `[n_pop × n_regions × n_yrs × n_seas × n_cat × n_sexes × n_fleets × n_sims]`.

- ISS:

  Integer sample sizes
  `[n_regions × n_yrs × n_seas × n_sexes × n_fleets × n_sims]`, read
  when `pop_specific = FALSE`.

- AgeingError:

  Ageing error matrices `[n_yrs × n_obs_ages × n_ages × n_sims]`,
  ignored when `age_or_len = 1`.

- comp_like:

  Integer vector `[n_fleets]` of the likelihood per fleet: `0`
  multinomial, `1` Dirichlet-multinomial, `2`-`4` the logistic-normal
  forms.

- ln_theta:

  Log overdispersion `[n_regions × n_sexes × n_fleets]`, read when
  `pop_specific = FALSE`.

- corr_pars:

  Logistic-normal correlation parameters
  `[n_regions × n_sexes × n_fleets × n_corr_pars]`.

- ln_theta_agg, corr_pars_agg:

  Their counterparts for the aggregated compositions, each of length
  `n_fleets`.

- comp_type:

  Integer matrix `[n_yrs × n_fleets]`: `0` aggregated across regions,
  `1` split by sex, `2` joint across sexes, `999` no data.

- n_sexes, n_pop, n_regions, n_cat:

  Model dimensions, `n_cat` being the number of ages or lengths.

- Obs:

  Observed composition container, dimensioned like `Exp` and written in
  place.

- pop_specific:

  Logical. `TRUE` simulates each population separately from
  population-specific inputs.

- ISS_pop:

  Population-specific sample sizes
  `[n_pop × n_regions × n_yrs × n_seas × n_sexes × n_fleets × n_sims]`,
  read when `pop_specific = TRUE`.

- pop_comp_like, pop_comp_type:

  The likelihood and aggregation structure for the population-specific
  compositions, shaped as their aggregate counterparts.

- ln_pop_theta:

  Log overdispersion `[n_pop × n_regions × n_sexes × n_fleets]`.

- pop_corr_pars:

  Logistic-normal correlation parameters
  `[n_pop × n_regions × n_sexes × n_fleets × n_corr_pars]`.

- ln_pop_theta_agg, pop_corr_pars_agg:

  Their counterparts for the population-specific aggregated
  compositions, each `[n_pop × n_fleets]`.

- age_or_len:

  Integer. `0` for age compositions, which take ageing error, `1` for
  length compositions, which do not.

## Value

`Obs` with the draws filled in at `[r, y, seas, , , f, sim]`, or
`[p, r, y, seas, , , f, sim]` under `pop_specific = TRUE`. Every other
slice is unchanged.

## Details

Joint compositions (`comp_type = 2`) apply ageing error across the
combined age by sex vector through the Kronecker product `diag(n_sexes)`
and `AgeingError`. Aggregated compositions (`comp_type = 0`) are drawn
only on the final region pass, from expected proportions marginalized
over regions and sexes. Under `pop_specific = TRUE` each population is
drawn separately from its own sample sizes, dispersion and correlations,
and that aggregation happens within a population.
