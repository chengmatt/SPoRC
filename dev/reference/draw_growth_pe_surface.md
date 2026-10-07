# One draw of a deviation surface from the process the estimation model penalizes

The reverse of
[`Get_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_PE_loglik.md).
A shared level is drawn once and written wherever it appears, at the sd
of the slot the penalty reads it at. Under iid every level is its own
normal. Under the random walk each level steps from the year before. A
first year the estimation model gives a diffuse start (`rw_init_sigma`)
keeps the replicate's own value, since that start leaves the level to
the data, as the recruitment and F walks do; under `NA` it starts at the
walk's own sd, as the penalty does. The 3D GMRF and the separable AR1
draw each population, region and sex's whole surface over years and
`bins` from the form's precision, conditional on the cells the map fixes
at zero, which is the density the penalty evaluates at them.

## Usage

``` r
draw_growth_pe_surface(
  PE_model,
  map,
  pe_pars,
  bins,
  fit_devs = NULL,
  n_cond = 0,
  rw_init_sigma = NA
)
```

## Arguments

- PE_model:

  Integer process error code, as `Get_PE_loglik` reads it.

- map:

  Integer array `[pop, region, year, bin, sex]` of the levels, `NA`
  where a cell is fixed.

- pe_pars:

  Array `[pop, region, slot, sex]`, the half of `growth_pe_pars` the
  form reads.

- bins:

  Bins the correlated forms run over.

- fit_devs:

  The fit's deviations, shaped as `map`, read in the first `n_cond`
  years only.

- n_cond:

  Integer. Years that keep the fit's deviations. Default `0` draws every
  year.

- rw_init_sigma:

  The sd the estimation model gives a walk's first year, or `NA`
  (default) for the walk's own.

## Value

Array shaped as `map`. A cell the map fixes keeps the fit's value, as
the penalty reads it, zero when no fit is given.

## Details

The first `n_cond` years keep the fit's deviations, and the years after
them are drawn given those: a random walk steps on from the fit's last
value, and the correlated forms draw from the precision conditional on
the fit's cells as well as the fixed ones.
