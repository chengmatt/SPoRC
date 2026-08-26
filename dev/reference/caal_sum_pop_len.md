# Sum a conditional age-at-length array across populations, for one length bin

As
[`caal_sum_pop`](https://chengmatt.github.io/SPoRC/dev/reference/caal_sum_pop.md),
for a single length bin.

## Usage

``` r
caal_sum_pop_len(arr, y, seas, l, f, n_pop, n_regions, n_ages, n_sexes)
```

## Arguments

- arr:

  Array indexed population, region, year, season, length, age, sex,
  fleet.

- y, seas, l, f:

  Year, season, length bin and fleet to extract.

- n_pop, n_regions, n_ages, n_sexes:

  Model dimensions.

## Value

An array indexed region, age, sex.
