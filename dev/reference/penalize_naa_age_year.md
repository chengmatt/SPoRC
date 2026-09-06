# Compare one region's age-by-year innovation surface

The age and year half of
[`Get_NAA_state_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_NAA_state_penalty.md),
split out so the region correlation can whiten its dim and then reuse
this unchanged for every structure, including the three-dimensional
field whose cohort term makes it non-separable.

## Usage

``` r
penalize_naa_age_year(eps_ya, sd_prs, NAA_re, pe, ny, na)
```

## Arguments

- eps_ya:

  Matrix `[year, age]` of innovations, already whitened across regions
  when a region correlation is active.

- sd_prs:

  Standard deviation for this population, region and sex.

- NAA_re:

  Integer structure code, as in
  [`Get_NAA_state_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_NAA_state_penalty.md).

- pe:

  Numeric vector of three correlation parameters on the unconstrained
  scale, read as age, year and cohort.

- ny, na:

  Number of active years and ages.

## Value

Scalar negative log likelihood.
