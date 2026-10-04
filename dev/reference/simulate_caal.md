# Simulate conditional age-at-length observations

Draws one age composition per length bin from the joint numbers at
length and age the fit builds for the same fleet: the catch or index at
each age spread over length by the size-age key, or under selectivity at
length the fish available at each age spread over length and selected
length by length. The row for a region, length bin and sex is summed
over populations and read through the fleet's ageing error onto the
observed ages, and the draw for bin \\l\\ is a multinomial (or
Dirichlet-multinomial) of `ISS[l]` fish across observed ages with that
row as the probability, which is the conditional \\P(a \mid l)\\ the fit
evaluates.

## Usage

``` r
simulate_caal(
  r,
  y,
  f,
  seas,
  sim,
  Joint,
  ISS,
  AgeingError,
  comp_like,
  ln_theta,
  ln_theta_agg,
  comp_type,
  n_sexes,
  n_regions,
  n_lens,
  Obs,
  CAAL_LenBinMap = NULL
)
```

## Arguments

- r, y, f, seas, sim:

  Region, year, fleet, season and replicate indices.

- Joint:

  Array `[pop, region, len, age, sex]` of the fleet's numbers at length
  and age in this year, season and replicate.

- ISS:

  Array `[region, year, season, length row, sex, fleet, sim]` of fish
  aged per length row. A zero skips the row.

- AgeingError:

  Array `[year, model_age, obs_age, sim]`.

- comp_like:

  Likelihood code per fleet (0 multinomial, 1 DM, 999 none).

- ln_theta:

  Array `[region, sex, fleet]` of DM log overdispersion.

- ln_theta_agg:

  Vector of aggregated DM log overdispersion per fleet.

- comp_type:

  Matrix `[year, fleet]` of composition type codes.

- n_sexes, n_regions, n_lens:

  Dimension sizes.

- Obs:

  Array `[region, year, season, length row, obs_age, sex, fleet, sim]`
  the draws are written into.

- CAAL_LenBinMap:

  Optional 0/1 matrix `[n_lens x n_caal_lens]` of the model length bins
  each length row covers. `NULL` (default) gives one row per model bin.

## Value

The updated `Obs` array.

## Details

Composition types follow `simulate_comps`: split by region and sex (1)
draws each sex separately, joint by sex (2) draws one sample across the
age by sex stack, and aggregated (0) pools regions and sexes and is
drawn once when the last region is reached. Only the multinomial and
Dirichlet multinomial families exist for CAAL.
