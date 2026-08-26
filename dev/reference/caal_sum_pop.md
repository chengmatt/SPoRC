# Sum a conditional age-at-length array across populations

The conditional age-at-length arrays carry a population dimension the
likelihood does not use, so it is summed away before the comparison. A
single population needs only a reshape, which avoids an apply over a
degenerate margin.

## Usage

``` r
caal_sum_pop(arr, y, seas, f, n_pop, n_regions, n_lens, n_ages, n_sexes)
```

## Arguments

- arr:

  Array indexed population, region, year, season, length, age, sex,
  fleet.

- y, seas, f:

  Year, season and fleet to extract.

- n_pop, n_regions, n_lens, n_ages, n_sexes:

  Model dimensions.

## Value

An array indexed region, length, age, sex.
