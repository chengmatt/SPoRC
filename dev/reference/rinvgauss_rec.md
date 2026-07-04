# Sample recruitment from an inverse-Gaussian distribution

Generates `sims` random recruitment values from an inverse-Gaussian
distribution whose parameters are estimated from a historical
recruitment vector using the method of moments. The arithmetic mean
\\\bar{R}\_a\\ and harmonic mean \\\bar{R}\_h\\ of `recruitment` are
used to derive the shape parameter \\\delta = 1 / (\bar{R}\_a /
\bar{R}\_h - 1)\\ and scale parameter \\\beta = \bar{R}\_a\\. Random
draws are generated via the Michael–Schucany–Haas (1976)
acceptance-mixture algorithm: a squared standard normal variate \\\psi\\
yields two candidate roots \\\omega\\ and \\\zeta\\, and a uniform draw
selects between them with probability \\\beta / (\beta + \omega)\\.

## Usage

``` r
rinvgauss_rec(sims, recruitment)
```

## Arguments

- sims:

  Integer. Number of random recruitment values to generate.

- recruitment:

  Numeric vector of historical recruitment values (must be strictly
  positive). Used to estimate \\\bar{R}\_a\\ and \\\bar{R}\_h\\.

## Value

Numeric vector of length `sims` of positive random variates drawn from
the fitted inverse-Gaussian distribution.

## Details

The inverse-Gaussian is appropriate for recruitment time series that are
right-skewed and strictly positive, and is used in SPoRC's resampling
recruitment option (`recruitment_opt = 999`) during projection years.
