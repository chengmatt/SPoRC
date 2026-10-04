# Sum a conditional age-at-length array across populations, for one length row

As
[`caal_sum_pop`](https://chengmatt.github.io/SPoRC/dev/reference/caal_sum_pop.md),
for a single length row.

## Usage

``` r
caal_sum_pop_len(
  arr,
  y,
  seas,
  l,
  f,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  caal_len_bin_map = NULL
)
```

## Arguments

- arr:

  Array indexed population, region, year, season, length, age, sex,
  fleet.

- y, seas, l, f:

  Year, season, length row and fleet to extract.

- n_pop, n_regions, n_ages, n_sexes:

  Model dimensions.

- caal_len_bin_map:

  Optional 0/1 matrix `[n_lens x n_caal_lens]`, the model length bins
  each age-at-length row covers. `NULL` (default) makes `l` a model
  length bin.

## Value

An array indexed region, age, sex.
