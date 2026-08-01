# Map fishing mortality parameters

Constructs the `ln_F_devs` and `ln_F_mean` factor maps, assigning unique
estimation indices to cells where catch data are used (`UseCatch == 1`)
and mapping cells without catch data to `NA`. This ensures that fishing
mortality parameters are only estimated for dimensions with observed
catch. `ln_F_devs` is resolved per region-year-season-fleet cell, while
`ln_F_mean` is resolved per region-season-fleet cell and is estimated
whenever that cell is fished in at least one year.

## Usage

``` r
do_Fmort_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$UseCatch` and `$data$UseCatch_pop` to be populated by
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

## Value

The input `input_list` with `$map$ln_F_devs` and `$map$ln_F_mean` set to
factor vectors. Cells with catch are assigned sequential integer
indices; cells without catch are `NA`.
