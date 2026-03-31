# Helper function used to check simulation object dimensions to ensure they are correct

Helper function used to check simulation object dimensions to ensure
they are correct

## Usage

``` r
check_sim_dimensions(
  x,
  n_regions = NULL,
  n_years = NULL,
  n_ages = NULL,
  n_lens = NULL,
  n_sexes = NULL,
  n_fish_fleets = NULL,
  n_srv_fleets = NULL,
  n_sims = NULL,
  what
)
```

## Arguments

- x:

  Object to evaluate

- n_regions:

  Number of regions

- n_years:

  Number of years

- n_ages:

  Number of ages

- n_lens:

  Number of lengths

- n_sexes:

  Number of sexes

- n_fish_fleets:

  Number of fishery fleets

- n_srv_fleets:

  Number of survey fleets

- n_sims:

  Number of simulations

- what:

  character specifying what to be evaluated
