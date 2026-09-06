# Expand a fleet-specific ageing error specification to its full array

Ageing error is a property of the sampling program, so a fishery that
ages its catch from otoliths and a survey that reads scales do not
misclassify the same way. `AgeingError_fish` and `AgeingError_srv` let
each fleet have its own matrix, defaulting to the shared `AgeingError`
so a model written before they existed behaves exactly as it did.

## Usage

``` r
expand_fleet_ageing_error(x, shared, n_fleets, what)
```

## Arguments

- x:

  `NULL`, a 3D array `[n_ages x n_obs_ages x n_fleets]` for a
  time-invariant fleet-specific matrix, or a 4D array
  `[n_years x n_ages x n_obs_ages x n_fleets]` for a time-varying one.

- shared:

  The shared `[n_years x n_ages x n_obs_ages]` array to fall back on,
  already expanded over years.

- n_fleets:

  Integer. Number of fleets.

- what:

  Character. Argument name, used in messages and errors.

## Value

Array `[n_years x n_ages x n_obs_ages x n_fleets]`.

## Details

Every fleet must land on the same observed age bins, because the
observed composition arrays have a single age dimension shared across
fleets.
