# Validate an internal-OSA `method`

The internal OSA path deliberately restricts
[`RTMB::oneStepPredict()`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)'s
`method` to the generic/Gaussian family. In particular the `"cdf"`
method is disallowed: it is numerically fragile for the discrete
(multinomial / count) likelihoods used here and can silently return
mis-calibrated residuals, so only `"oneStepGeneric"`,
`"oneStepGaussian"`, and `"oneStepGaussianOffMode"` are accepted.

## Usage

``` r
validate_osa_method(method)
```

## Arguments

- method:

  Character scalar method name.

## Value

`method` invisibly, if valid; otherwise an error is raised.
