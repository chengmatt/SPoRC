# How many observed ages the at-age data sources are recorded on

At-age observations are read through the ageing error matrix, so they
sit on its columns, the same observed ages the age compositions use. A
list with no ageing error yet, as the mapping unit tests assemble, has
them on the model ages.

## Usage

``` r
at_age_n_obs_ages(input_list)
```

## Arguments

- input_list:

  Named list with `$data`.

## Value

An integer.
