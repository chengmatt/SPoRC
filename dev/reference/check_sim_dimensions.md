# Validate Dimensions of Simulation Input Objects

Internal utility that verifies the dimensional structure of simulation
input objects used in the operating model. The function checks that the
supplied object `x` has the expected dimensions for the object type
specified by `what`. If the dimensions do not match the required
structure, the function stops with an informative error.

## Usage

``` r
check_sim_dimensions(
  x,
  n_pop = NULL,
  n_regions = NULL,
  n_years = NULL,
  n_ages = NULL,
  n_lens = NULL,
  n_sexes = NULL,
  n_seas = NULL,
  n_fish_fleets = NULL,
  n_srv_fleets = NULL,
  n_sims = NULL,
  what
)
```

## Arguments

- x:

  Object to evaluate. Typically a numeric array (or vector) whose
  dimensions must match the expected structure associated with `what`.

- n_pop:

  Integer. Number of populations.

- n_regions:

  Integer. Number of spatial regions.

- n_years:

  Integer. Number of model years.

- n_ages:

  Integer. Number of age classes.

- n_lens:

  Integer. Number of length bins.

- n_sexes:

  Integer. Number of sexes.

- n_seas:

  Integer. Number of seasons.

- n_fish_fleets:

  Integer. Number of fishery fleets.

- n_srv_fleets:

  Integer. Number of survey fleets.

- n_sims:

  Integer. Number of stochastic simulations.

- what:

  Character string identifying the type of object to validate. This
  determines the expected dimension ordering and structure.

## Value

No return value. The function is used for validation and terminates with
an error if the dimensions of `x` do not match the expected structure.

## Details

Simulation inputs generally include an additional dimension representing
the number of stochastic simulations (`n_sims`). The position of this
dimension depends on the object type but is typically the final
dimension of the array.

The function validates the dimensional structure of several classes of
simulation inputs used in the operating model:

- Biological processes (e.g., weight-at-age, maturity, natural
  mortality)

- Fishery processes (e.g., fishing mortality, selectivity, catchability)

- Survey processes (e.g., selectivity, catchability)

- Observation models for indices and composition data

- Recruitment and demographic processes

- Conventional tagging parameters

Arrays generally follow the dimension ordering


    population × region × year × season × age × sex × fleet × simulation

although only the dimensions relevant to a given object are included.
