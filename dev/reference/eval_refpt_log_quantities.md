# Evaluate reference point quantities at a trial parameter vector

Evaluate reference point quantities at a trial parameter vector

## Usage

``` r
eval_refpt_log_quantities(
  obj,
  p,
  refpt_args,
  extra_quantities = NULL,
  keep = NULL
)
```

## Arguments

- obj:

  Fitted RTMB model object.

- p:

  Numeric vector ordered as `obj$env$last.par.best`.

- refpt_args:

  List forwarded to
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md).

- extra_quantities:

  Optional `function(rep, refpts)` returning a named vector of extra
  positive quantities.

- keep:

  Optional character vector of quantity names to retain.

## Value

Named numeric vector on the log scale.
