# Extract the (AD) values slot from an OSA observation

Returns the observation values with their AD class preserved, so the
single bin being peeled by `oneStepPredict` remains differentiable.
Contrast with
[`osa_extract_values`](https://chengmatt.github.io/SPoRC/dev/reference/osa_extract_values.md),
which deliberately detaches the AD class to freeze the running
remainder.

## Usage

``` r
osa_extract_x(xobs)
```

## Arguments

- xobs:

  Either an object of class `"osa"` or a plain numeric vector.

## Value

The AD (or numeric) values of the observation.
