# Combine a parameter function and a data list for RTMB

Returns a closure that calls `f(p, d)`, allowing the data list to be
fixed at construction time so that `RTMB::MakeADFun` receives a
single-argument objective function of the form `function(p)`.

## Usage

``` r
cmb(f, d)
```

## Arguments

- f:

  Function with signature `function(pars, data)`, typically
  [`SPoRC_rtmb`](https://chengmatt.github.io/SPoRC/dev/reference/SPoRC_rtmb.md).

- d:

  Named list of model data passed as the second argument to `f` on every
  call.

## Value

A single-argument function `function(p)` equivalent to `f(p, d)`.
