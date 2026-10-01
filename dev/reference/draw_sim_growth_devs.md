# Draw growth deviations for every replicate and rebuild growth

Each varying growth parameter's series and the semi-parametric surface
are drawn by
[`draw_growth_pe_surface`](https://chengmatt.github.io/SPoRC/dev/reference/draw_growth_pe_surface.md)
from the form the estimation model penalizes them under, at the fitted
process error parameters: a parameter's sd in the time-varying half of
`growth_pe_pars` and the surface's in the semi-parametric half, as
`Get_PE_loglik` reads them. The fit's maps decide what is drawn: a cell
the map fixes stays at zero, a shared level takes one draw, and years
past the fit are active. Growth is then rebuilt through
`derive_sim_growth`; under cohort growth that builds only the years
before the propagation starts, and the annual cycle advances the rest
from each replicate's own numbers at age.

## Usage

``` r
draw_sim_growth_devs(sim_env)
```

## Arguments

- sim_env:

  Simulation environment holding what `Setup_Sim_Growth_RE` stored,
  `n_yrs` and `n_sims`.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
