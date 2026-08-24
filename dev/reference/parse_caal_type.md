# Translate a CAAL composition type specification into a year by fleet matrix

Uses the same `"CompType_Year_x-y_Fleet_z"` string convention as the
marginal compositions, so a CAAL type is specified exactly the way an
age or length composition type is. `"terminal"` in place of the upper
year runs the block to the last model year.

## Usage

``` r
parse_caal_type(CAAL_Type, n_yrs, n_fleets, what)
```

## Arguments

- CAAL_Type:

  Character vector of composition type specifications.

- n_yrs:

  Number of model years.

- n_fleets:

  Number of fleets.

- what:

  Name used in error messages.

## Value

Integer matrix `[year x fleet]` of composition type codes (0 aggregated,
1 split by region and sex, 2 joint by sex and split by region, 999
none).
