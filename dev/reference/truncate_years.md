# Truncate an array's year dimension

Keeps the first `n_years` of one dim without caring how many dims the
array has, so arrays that have gained a dim still work.

## Usage

``` r
truncate_years(arr, n_years, yr_dim = 3)
```

## Arguments

- arr:

  Array to truncate. `NULL` passes through.

- n_years:

  Number of years to keep.

- yr_dim:

  Position of the year dim.

## Value

`arr` with years truncated.
