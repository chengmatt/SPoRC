# Standard deviation for one at-age observation

An at-age observation may have its own reported standard error, an
estimated component, or both, matching what the aggregated index data
sources allow. The parameter alone is the default and is what a data
source with no reported errors means.

## Usage

``` r
at_age_obs_sd(se, extra, form)
```

## Arguments

- se:

  Reported standard errors for the ages in one cell.

- extra:

  Estimated component for the same ages, on the natural scale.

- form:

  Integer. `0` the parameter alone, `1` the reported errors alone, `2`
  additive, `3` in quadrature.

## Value

A vector of standard deviations the length of `extra`.
