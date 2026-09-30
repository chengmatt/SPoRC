# Parameters each replicate's operating model runs on

Conditional gives every replicate the fitted parameters, so one truth
covers them all and only the data change. Joint gives each replicate its
own draw, so each has its own truth.

## Usage

``` r
sim_draw_views(
  sim_type,
  n_sims,
  fit_rep,
  parameters,
  mapping,
  sd_rep,
  random,
  obj
)
```

## Arguments

- sim_type:

  Either `"conditional"` or `"joint"`.

- n_sims:

  Number of replicates.

- fit_rep:

  Report list from the fitted model.

- parameters:

  Parameter list the model was built with.

- mapping:

  Factor maps the model was built with.

- sd_rep:

  `sdreport` object from the fitted model.

- random:

  Character vector of random effect names.

- obj:

  Fitted object, needed only under `"joint"`.

## Value

`reps` and `pars`, one per replicate, plus `fit_pars` at the fit.

## Details

The draw is taken at the fitted parameter vector, the order the joint
precision is in, and `parList` puts the map back. A fit with no random
effects has no joint precision, so the fixed effect covariance is
inverted in its place.
