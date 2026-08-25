# Gather the arguments [`compute_mortality_year`](https://chengmatt.github.io/SPoRC/dev/reference/compute_mortality_year.md) needs from the model frame

Gather the arguments
[`compute_mortality_year`](https://chengmatt.github.io/SPoRC/dev/reference/compute_mortality_year.md)
needs from the model frame

## Usage

``` r
mortality_args_from_model(env = parent.frame())
```

## Arguments

- env:

  Environment holding the unpacked data and parameters, i.e. the
  `SPoRC_rtmb` frame after
  [`RTMB::getAll`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).
  Defaults to the caller.

## Value

A named list of arguments for
[`compute_mortality_year`](https://chengmatt.github.io/SPoRC/dev/reference/compute_mortality_year.md).
