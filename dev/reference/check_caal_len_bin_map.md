# Validate the length rows of the conditional age-at-length data

Each column is one length row of the age-at-length data and holds a 1
for every model length bin that row covers, so a row can span several
bins.

## Usage

``` r
check_caal_len_bin_map(x, n_model_bins, what)
```

## Arguments

- x:

  The matrix to check.

- n_model_bins:

  Integer. Number of model length bins, the required row count.

- what:

  Character. Argument name, used in messages.

## Value

`x`, as a matrix.
