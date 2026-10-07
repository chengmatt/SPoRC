# A selectivity deviation map laid out by region and fleet

The units
[`Get_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_PE_loglik.md)
divides a shared deviation's penalty between, one row per region and
fleet whose penalty is evaluated. A deviation copied to another fleet
under `"est_shared_f_x"` is then one parameter held by both fleets and
penalized once, as one shared over regions already is.

## Usage

``` r
sel_devs_by_unit(map, cont_tv, pe_wt)
```

## Arguments

- map:

  Array `[region, year, bin, sex, fleet]` of deviation levels.

- cont_tv:

  Integer array `[region, fleet]` of time variation codes.

- pe_wt:

  Penalty weight by fleet; a fleet at zero is not penalized.

## Value

Array `[n_regions * n_fleets, year, bin, sex, 1]`, unit
`r + (f - 1) * n_regions`, `NA` where a unit is not penalized.
