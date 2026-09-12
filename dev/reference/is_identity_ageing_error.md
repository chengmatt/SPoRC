# Is one fleet's ageing error the identity in every year?

An identity matrix reads every age as itself, so a fleet with one is
predicted on the model's own ages without passing through the map.

## Usage

``` r
is_identity_ageing_error(ageing_error)
```

## Arguments

- ageing_error:

  Array `[n_years, n_ages, n_obs_ages]` for one fleet.

## Value

`TRUE` when every year is the identity matrix.
