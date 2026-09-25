# Starting mean of a covariate on its link scale

The observed mean through the link: itself under identity, its log under
the log link, its logit or complementary log-log under those (the mean
kept inside 0.02 to 0.98 first). A lognormal takes the mean log instead.
Zero when nothing is observed or the value is not finite.

## Usage

``` r
dsem_cov_link_mean(y, family, link)
```

## Arguments

- y:

  Observations with NA for missing years.

- family:

  Family code.

- link:

  Link code.

## Value

Scalar.
