# Draw one replicate's selectivity deviations past the fit and rebuild those years

The deviations are drawn on from the replicate's fitted ones under the
fitted process and rebuilt with `Get_Selex_Array` over the fitted and
projection years, the projection years reading the terminal year's
blocks and model as the estimation model's projection does. A form
centered over all its years would put the projection years in the mean
and shift every fitted year, so each region, sex and fleet's projection
years are scaled by what takes the longer build's terminal fitted year
back to the fit's.

## Usage

``` r
draw_sel_devs_past_fit(sim_env, type, data, pars, args, selex, n_proj_om, sim)
```

## Arguments

- sim_env:

  Simulation environment.

- type:

  `"fish"`, `"ret"` or `"srv"`.

- data, pars:

  The fit's data list and this replicate's parameters.

- args:

  The `Get_Selex_Array` arguments the fitted years were rebuilt with,
  their deviations the replicate's.

- selex:

  That rebuild's result.

- n_proj_om:

  Years the operating model runs past the fit.

- sim:

  Replicate.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
