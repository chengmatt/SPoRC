# Hold a natural mortality array across seasons

M is an instantaneous rate per year on a
`[n_pop x n_regions x n_years x n_seas x n_ages x n_sexes]` grid. An
array without the season dim gets kept at one rate across seasons, which
is the model it came from. An array that already has it passes through.

## Usage

``` r
expand_natmort_seasons(x, n_seas, seas_dim = 4, n_dim = 6)
```

## Arguments

- x:

  Natural mortality array, with or without seasons. `NULL` passes
  through.

- n_seas:

  Number of seasons.

- seas_dim:

  Where the season dim sits. 4 is the model's layout; the recruitment
  routines use a smaller grid with seasons third.

- n_dim:

  Number of dims once seasons are there.

## Value

An `n_dim` array with seasons in `seas_dim`.
