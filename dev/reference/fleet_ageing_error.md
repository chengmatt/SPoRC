# Read a fitted model's fleet-specific ageing error, falling back on the shared one

The post-fit diagnostics have to reproduce the expected compositions the
likelihood built, which means reading the same ageing error the
objective read. Models fitted before `AgeingError_fish` and
`AgeingError_srv` existed carry only the shared matrix, so it is
replicated across the fleets and the diagnostics come out exactly as
they did.

## Usage

``` r
fleet_ageing_error(data, shared, which)
```

## Arguments

- data:

  Data list from the fitted model.

- shared:

  The shared `[n_years x n_ages x n_obs_ages]` array.

- which:

  Either `"fish"` or `"srv"`.

## Value

Array `[n_years x n_ages x n_obs_ages x n_fleets]`.
