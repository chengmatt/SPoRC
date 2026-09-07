# Process error density for one series of recruitment deviations

The recruitment deviations of one population and region are either
independent draws about a supplied mean, a random walk, or an AR1
process. The walk and the AR1 always step from the previous year,
whether or not that year's deviation is estimated: every year has a
recruitment, so a deviation fixed at a value is a year the walk passes
through rather than a gap in the series. A year whose own deviation is
fixed contributes no density.

## Usage

``` r
get_recdev_pe_nLL(devs, is_est, sigma, dev_mu, PE_model, rho = 0, init_sd = 5)
```

## Arguments

- devs:

  Numeric vector of deviations for one population and region, by year.

- is_est:

  Numeric vector the same length, `1` where the deviation is estimated
  and `0` where it is fixed.

- sigma:

  Numeric vector the same length of standard deviations, so the early
  and late regimes can differ. A step reads the standard deviation of
  the year it lands on.

- dev_mu:

  Numeric vector the same length of prior means. Only read when
  `PE_model = 1`, since a walk's mean is the previous deviation.

- PE_model:

  Integer. `1` independent, `2` random walk, `3` AR1.

- rho:

  AR1 correlation, already on the natural scale. Only read when
  `PE_model = 3`.

- init_sd:

  Standard deviation given to year one of a random walk. Default `5`,
  which leaves the level effectively free. `NA` instead starts the walk
  at zero under its own sigma.

## Value

Numeric vector of negative log-likelihood contributions, zero where the
deviation is fixed.

## Details

A walk has no stationary distribution to start from, so year one is
given a diffuse normal. An AR1 starts from its stationary marginal
standard deviation \\\sigma / \sqrt{1 - \rho^2}\\. Fixing year one
instead leaves the series with no penalty on its level at all, which is
what SAM's flat prior on the first year amounts to.
