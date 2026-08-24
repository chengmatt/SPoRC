# Simulate conditional age-at-length observations

Draws one age composition per length bin from the joint distribution of
length and age implied by the size-age transition matrix and the true
numbers at age (catch at age for the fishery, index at age for the
survey). The joint for a region, length bin and sex is \\P(l \mid a)
N_a\\ summed over populations, and the draw for bin \\l\\ is a
multinomial (or Dirichlet-multinomial) of `ISS[l]` fish across ages with
that row as the probability, which is the conditional \\P(a \mid l)\\ by
construction. Ageing error is applied to the drawn counts the same way
[`simulate_comps`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_comps.md)
applies it to marginal age compositions.

## Usage

``` r
simulate_caal(
  r,
  y,
  f,
  seas,
  sim,
  SizeAgeTrans,
  AtAge,
  ISS,
  AgeingError,
  comp_like,
  ln_theta,
  ln_theta_agg,
  comp_type,
  n_sexes,
  n_regions,
  n_lens,
  Obs
)
```

## Arguments

- r, y, f, seas, sim:

  Region, year, fleet, season and replicate indices.

- SizeAgeTrans:

  Array `[pop, region, year, season, len, age, sex, sim]` of \\P(l \mid
  a)\\.

- AtAge:

  Array `[pop, region, year, season, age, sex, fleet, sim]` of true
  numbers at age for this fleet type.

- ISS:

  Array `[region, year, season, len, sex, fleet, sim]` of fish aged per
  length bin. A zero skips the bin.

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

  Array `[region, year, season, len, obs_age, sex, fleet, sim]` the
  draws are written into.

## Value

The updated `Obs` array.

## Details

Composition types follow `simulate_comps`: split by region and sex (1)
draws each sex separately, joint by sex (2) draws one sample across the
age by sex stack, and aggregated (0) pools regions and sexes and is
drawn once when the last region is reached. Only the multinomial and
Dirichlet multinomial families exist for CAAL.
