# Extract the lower/upper CDF indicators from an OSA observation

Returns the `cdf_lower` / `cdf_upper` data indicators that the `"cdf"`
method of
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
toggles. When `xobs` is a plain numeric vector (ordinary fitting) these
are all zero, so the CDF terms contribute nothing and the density
reduces to the ordinary likelihood.

## Usage

``` r
osa_extract_cdf(xobs, n)
```

## Arguments

- xobs:

  An `"osa"` object or plain numeric vector.

- n:

  Length for the default all-zero indicators.

## Value

A list with numeric/AD vectors `lower` and `upper`.
