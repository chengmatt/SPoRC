# One draw of a fleet type's selectivity deviations

Each region and fleet's surface is drawn by
[`draw_growth_pe_surface`](https://chengmatt.github.io/SPoRC/dev/reference/draw_growth_pe_surface.md)
under that cell's process, the region in the place it reads as
population. A level shared over regions, sexes, bins or fleets then
takes the value drawn at its first cell, where the penalty reads its sd,
which also fills the bins the correlated forms leave out of their shared
groups. A level shared over units whose sds differ is drawn at the first
unit's.

## Usage

``` r
draw_sel_dev_surface(
  PE_codes,
  map,
  pe_pars,
  fit_devs,
  bins,
  n_cond,
  pe_wt = NULL,
  rw_init_sigma = NULL
)
```

## Arguments

- PE_codes:

  Integer array `[n_regions, n_fleets]` of process codes, or a vector
  `[n_fleets]` shared by every region.

- map:

  Array `[n_regions, n_yrs, n_bins, n_sexes, n_fleets]` of levels, `NA`
  where fixed.

- pe_pars:

  Array `[n_regions, slot, n_sexes, n_fleets]`.

- fit_devs:

  The fit's deviations, shaped as `map`.

- bins:

  Bins the correlated forms run over.

- n_cond:

  Leading years that keep the fit's deviations.

- pe_wt:

  Penalty weight by fleet; a fleet at zero has no penalty and keeps the
  fit's.

- rw_init_sigma:

  The sd the estimation model gives a walk's first year, by fleet, `NA`
  for the walk's own. `NULL` (default) is `NA` for every fleet.

## Value

Array shaped as `map`.
