# Extract frozen numeric values from an OSA observation

Returns the numeric *values* of an observation slice, detached from the
AD tape. This is the R analogue of WHAM's `asDouble()`: it is used to
build the running composition remainder so that, when
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
peels and perturbs a single bin, the "all other bins" count does not
move with the perturbation. Without this the running remainder can be
driven negative on candidate values, producing `NaN` density evaluations
(most visibly with random effects switched off, where the generic
integrator sweeps the full support).

## Usage

``` r
osa_extract_values(xobs)
```

## Arguments

- xobs:

  Either an object of class `"osa"` (with slots `@x` and `@keep`)
  supplied by `oneStepPredict`, or a plain numeric vector used during
  ordinary fitting.

## Value

A plain numeric vector of observed counts.
