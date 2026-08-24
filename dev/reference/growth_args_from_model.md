# Assemble the growth module's arguments from the model's data and parameters

The growth settings are the same for every call, whether the whole
series is built up front or one year at a time from inside the
population loop, so they are gathered once here rather than written out
at each call site.

## Usage

``` r
growth_args_from_model(env = parent.frame())
```

## Arguments

- env:

  Environment holding the unpacked data and parameters, i.e. the
  `SPoRC_rtmb` frame after
  [`RTMB::getAll`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).
  Defaults to the caller.

## Value

A named list of arguments for
[`Get_Growth`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth.md)
and
[`Get_Growth_Year`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth_Year.md).
