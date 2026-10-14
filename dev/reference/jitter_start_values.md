# Draw jittered starting values for one jitter iteration

[`nlminb()`](https://rdrr.io/r/stats/nlminb.html) only ever sees the
fixed effects, so a jittered random effect has to reach the model
another way. TMB starts the inner Laplace solve from `obj$env$last.par`,
so the draws are returned separately and written there by the caller.

## Usage

``` r
jitter_start_values(obj, par_vec, sd, jitter_random)
```

## Arguments

- obj:

  An RTMB model object from `MakeADFun`.

- par_vec:

  Starting values, either the fixed-effect vector (`length(obj$par)`) or
  the joint fixed and random vector (`length(obj$env$par)`), or `NULL`
  for the model's own start.

- sd:

  Standard deviation of the additive normal draws.

- jitter_random:

  Whether the random effects are perturbed as well.

## Value

A list with `fixed`, the jittered fixed-effect vector; `random`, the
random-effect starting values in `obj$env$random` order, jittered when
`jitter_random = TRUE` and taken straight from `par_vec` otherwise, and
of length zero when the model has no random effects; and `joint`, the
two combined into one vector in `obj$env$par` order.
