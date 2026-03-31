# Simulate Age or Length Compositions

Generates Observed fish compositions by age or length for a given
region, year, fleet, and simulation iteration. Supports multinomial,
Dirichlet-multinomial, and logistic-normal likelihoods, with optional
ageing error applied.

## Usage

``` r
simulate_comps(
  r,
  y,
  f,
  sim,
  Exp,
  ISS,
  AgeingError,
  comp_like,
  ln_theta,
  corr_pars,
  ln_theta_agg,
  corr_pars_agg,
  comp_type,
  n_sexes,
  n_regions,
  n_cat,
  Obs,
  age_or_len = 0
)
```

## Arguments

- r:

  Integer. Region index.

- y:

  Integer. Year index.

- f:

  Integer. Fleet index.

- sim:

  Integer. Simulation iteration index.

- Exp:

  Array. Expected compositions (age or length) with dimensions \[region,
  year, category, sex, fleet, sim\].

- ISS:

  Array. Sample size (integer) for the Observed compositions with
  dimensions \[region, year, sex, fleet, sim\].

- AgeingError:

  Array. Ageing error matrix for each year, dimensions \[year, category,
  category, sim\].

- comp_like:

  Integer vector. Composition likelihood type per fleet: 0 =
  multinomial, 1 = Dirichlet-multinomial, 2-4 = logistic-normal.

- ln_theta:

  Array. Log-variance parameter for compositions per region, sex, and
  fleet, dimensions \[region, sex, fleet\].

- corr_pars:

  Array. Correlation parameters for logistic-normal likelihood,
  dimensions \[region, sex, fleet, ?\].

- ln_theta_agg:

  Numeric vector. Log-variance parameter for aggregated compositions per
  fleet.

- corr_pars_agg:

  Numeric vector. Correlation parameters for aggregated logistic-normal
  likelihood per fleet.

- comp_type:

  Integer array. Composition type: 0 = aggregated across regions, 1 =
  split by sex, 2 = joint across sexes, dimensions \[year, fleet\].

- n_sexes:

  Integer. Number of sexes.

- n_regions:

  Integer. Number of regions.

- n_cat:

  Integer. Number of categories (ages or lengths).

- Obs:

  Array. Observed compositions array to fill, same dimensions as
  \`Exp\`.

- age_or_len:

  Integer. Flag to indicate if ageing error should be applied: 0 = apply
  ageing error (for ages), 1 = do not apply (for lengths).

## Value

Array of Observed compositions with the same dimensions as \`Obs\`,
updated with simulated Observations.

## Details

The function handles three cases based on \`comp_type\`: 1. Split by sex
(comp_type = 1): compositions are simulated separately for each sex. 2.
Joint compositions across sexes (comp_type = 2): compositions simulated
jointly and multiplied by a kronecker matrix for logistic-normal or
Dirichlet-multinomial likelihoods. 3. Aggregated across regions
(comp_type = 0): only applied in the last region and averages across
regions and sexes.

The function normalizes expected compositions, applies the selected
likelihood (\`comp_like\`), and multiplies by \`AgeingError\` when
applicable.
