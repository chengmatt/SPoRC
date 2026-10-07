# Standardized at-age residuals for the fleets whose ages are correlated

A correlated fleet cannot draw each age on its own in the annual cycle,
so its standardized residuals, each age's log-scale residual over that
age's sd, are drawn here for every replicate under the form the
estimation model evaluates (`get_at_age_nLL`, `get_at_age_2dar1_nLL`).
Under `"1dar1"` and `"us"` the ages a year observes are drawn together,
at a correlation of `rho^|age gap|` or the unstructured matrix subset to
those ages. Under `"2dar1"` the region, season and sex's whole block of
observed years by observed ages is one draw, an AR1 over the sequence of
observed years by an AR1 over ages. The residuals do not depend on the
population, so drawing them before the first year changes nothing else.

## Usage

``` r
draw_sim_at_age_corr(sim_env)
```

## Arguments

- sim_env:

  Simulation environment holding the at-age use flags and what
  [`sim_at_age_corr_setup`](https://chengmatt.github.io/SPoRC/dev/reference/sim_at_age_corr_setup.md)
  stored.

## Value

`invisible(NULL)`; `sim_env` gains `<data source>_std_resid`, shaped
like the use flags with the replicate dim last, for each data source
with a correlated fleet.
