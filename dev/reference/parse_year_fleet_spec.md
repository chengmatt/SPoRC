# Parse a year-by-fleet specification into a matrix

Settings that change part way through a series are given as strings of
the form `Value_Year_x-y_Fleet_f`, where `x` and `y` are positions in
the year vector rather than calendar years and `y` may be `"terminal"`.
A vector of them describes a whole model, one entry per block, and every
year of every fleet must be covered by some entry.

## Usage

``` r
parse_year_fleet_spec(spec, arg_name, n_fleets, n_yrs, codes, check = NULL)
```

## Arguments

- spec:

  Character vector of specifications.

- arg_name:

  Name of the argument being parsed, used in error messages.

- n_fleets:

  Number of fleets the matrix must cover.

- n_yrs:

  Number of years the matrix must cover.

- codes:

  Named numeric vector mapping each accepted value to the code the model
  reads. The names are the accepted vocabulary, so a stream with its own
  set of values passes its own mapping.

- check:

  Optional function called as `check(value, fleet)` for each entry,
  returning a character string to stop with or `NULL` to accept. Used
  for constraints that depend on other settings, such as a likelihood
  that cannot take an aggregated observation.

## Value

Numeric matrix `[n_yrs, n_fleets]` of codes.

## Details

Later entries overwrite earlier ones where they overlap, so a general
setting can be given first and a period carved out of it afterwards.
