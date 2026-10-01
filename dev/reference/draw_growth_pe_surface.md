# One draw of a deviation surface from the process the estimation model penalizes

The reverse of
[`Get_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_PE_loglik.md).
A shared level is drawn once and written wherever it appears, at the sd
of the slot the penalty reads it at. Under iid every level is its own
normal. Under the random walk each level steps from the year before; the
first year starts at the walk's own sd, as the recruitment walk does in
the operating model, since the diffuse start the estimation model gives
it only leaves the level free. The 3D GMRF and the separable AR1 draw
each population, region and sex's whole surface over years and `bins`
from the form's precision, conditional on the cells the map fixes at
zero, which is the density the penalty evaluates at them.

## Usage

``` r
draw_growth_pe_surface(PE_model, map, pe_pars, bins)
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

## Value

Array shaped as `map`, zero where the map is `NA`.
