# Penalty on catchability deviations

Annual deviations from a fleet's block catchability, taken as
independent, a random walk or an ar1. A fleet whose deviations a dsem
has taken over reads no penalty here, since the mirror blanks those
cells.

## Usage

``` r
Get_q_dev_penalty(
  ln_q_devs,
  ln_sigma_q,
  q_rho,
  q_model,
  map_ln_q_devs = NULL,
  q_rw_init_sigma = NA
)
```

## Arguments

- ln_q_devs:

  Array `[n_regions, n_yrs_total, n_fleets]` of log-scale catchability
  deviations.

- ln_sigma_q:

  Array `[n_regions, n_fleets]` of log-scale deviation standard
  deviations.

- q_rho:

  Array `[n_regions, n_fleets]` of unconstrained ar1 correlations, read
  only where `q_model` is ar1.

- q_model:

  Integer vector `[n_fleets]`: 1 none, 2 iid, 3 random walk, 4 ar1, 5
  dsem.

- map_ln_q_devs:

  Numeric array of map levels the same shape as `ln_q_devs`, `NA` where
  a cell is not penalized.

- q_rw_init_sigma:

  Standard deviation of the first estimated year of a random walk, or
  `NA` to start it at zero under its own sigma.

## Value

Negative log density, summed over regions, years and fleets.
