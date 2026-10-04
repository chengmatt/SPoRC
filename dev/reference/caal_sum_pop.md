# Sum a conditional age-at-length array across populations

The conditional age-at-length arrays have a population dimension the
likelihood does not use, so it is summed away before the comparison. A
single population needs only a reshape, which avoids an apply over a
degenerate dim. With `caal_len_bin_map`, each age-at-length row then
sums the model length bins it covers.

## Usage

``` r
caal_sum_pop(
  arr,
  y,
  seas,
  f,
  n_pop,
  n_regions,
  n_lens,
  n_ages,
  n_sexes,
  caal_len_bin_map = NULL
)
```

## Arguments

- arr:

  Array indexed population, region, year, season, length, age, sex,
  fleet.

- y, seas, f:

  Year, season and fleet to extract.

- n_pop, n_regions, n_lens, n_ages, n_sexes:

  Model dimensions.

- caal_len_bin_map:

  Optional 0/1 matrix `[n_lens x n_caal_lens]`, the model length bins
  each age-at-length row covers. `NULL` (default) keeps one row per
  model bin.

## Value

An array indexed region, length row, age, sex.
