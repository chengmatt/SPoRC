# Dimensions of one at-age observation array

Every at-age data source is stored region by year by season by age by
sex by fleet, with a leading population dimension for the
population-specific form. This is the layout of the prediction arrays
the likelihood reads, so the two line up dim for dim.

## Usage

``` r
at_age_dims(input_list, fleet_field, pop = FALSE)
```

## Arguments

- input_list:

  Named list with `$data`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- pop:

  Logical. `TRUE` for the population-specific data source.

## Value

An integer vector of dimensions.
