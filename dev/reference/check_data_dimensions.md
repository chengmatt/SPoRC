# Validate Dimensions of Model Input Objects

Internal utility that verifies the dimensional structure of model input
objects. The function checks that the supplied object `x` has the
expected dimensions for the specified object type `what`. Expected
dimension ordering depends on the data structure used in the population
model (e.g., population, region, year, season, age, sex, fleet).

## Usage

``` r
check_data_dimensions(
  x,
  n_pop = NULL,
  n_regions = NULL,
  n_years = NULL,
  n_ages = NULL,
  n_seas = NULL,
  n_lens = NULL,
  n_sexes = NULL,
  n_fish_fleets = NULL,
  n_srv_fleets = NULL,
  conv_tag_max_liberty = NULL,
  n_conv_tag_cohorts = NULL,
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

- n_seas:

  Integer. Number of seasons.

- n_lens:

  Integer. Number of length bins.

- n_sexes:

  Integer. Number of sexes.

- n_fish_fleets:

  Integer. Number of fishery fleets.

- n_srv_fleets:

  Integer. Number of survey fleets.

- conv_tag_max_liberty:

  Integer. Maximum tag liberty for conventional tagging data (number of
  seasons between release and recapture).

- n_conv_tag_cohorts:

  Integer. Number of conventional tagging cohorts.

- what:

  Character string identifying the type of object being validated. This
  determines the expected dimension ordering (e.g., biological inputs,
  fishery observations, survey observations, or tagging data).

## Value

No return value. The function is used for validation and will terminate
execution with an error if the dimensions of `x` do not match the
expected structure.

## Details

The argument `what` determines which dimension template is applied. If
the dimensions of `x` do not match the required structure, the function
stops with an informative error message.

The function supports validation of several model input classes
including:

- Biological quantities (e.g., weight-at-age, maturity, natural
  mortality)

- Movement and spatial processes

- Fishery observations (catch, indices, composition data)

- Survey observations (indices and compositions)

- Conventional tagging data and recaptures

Some objects (e.g., ageing-error matrices or age compositions) may
contain additional dimensions representing observed ages or length bins.
These dimensions are not always validated explicitly because they depend
on the structure of the observed data.
