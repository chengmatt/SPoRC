# Joint covariance of the fitted parameters and states

Selectivity is read from the terminal year, so a fit with random effects
needs the joint precision rather than the fixed effect block. The
precision is left sparse and never inverted.

## Usage

``` r
refpt_par_cov(obj, sd_rep = NULL)
```

## Arguments

- obj:

  Fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- sd_rep:

  Optional `sdreport` to reuse. Needs `jointPrecision` when the fit has
  random effects.

## Value

List with `kind` (`"cov"` or `"prec"`), `mat`, and `n`.
