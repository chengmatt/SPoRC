# Expand an at-age aggregation setting to a year by fleet matrix

The at-age streams accept their aggregation either as a bare value,
standing for the whole series, or as year and fleet specifications in
the same form the composition streams take. Both arrive here and leave
as a matrix, so everything downstream reads one shape.

## Usage

``` r
at_age_type_matrix(type, n_fleets, n_yrs, arg_name = "at-age Type")
```

## Arguments

- type:

  Character. Bare values, one for all fleets or one per fleet, or
  `Value_Year_x-y_Fleet_f` specifications.

- n_fleets:

  Number of fleets.

- n_yrs:

  Number of years.

- arg_name:

  Name of the argument being set, used in error messages.

## Value

Numeric matrix `[n_yrs, n_fleets]` of codes.
