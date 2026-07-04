# Map fishing mortality deviation parameters

Constructs the `ln_F_devs` factor map, assigning unique estimation
indices to region–year–season–fleet cells where catch data are used
(`UseCatch == 1`) and mapping cells without catch data to `NA`. This
ensures that `ln_F_devs` parameters are only estimated for dimensions
with observed catch.

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

The input `input_list` with `$map$ln_F_devs` set to a factor vector.
Cells with catch are assigned sequential integer indices; cells without
catch are `NA`.
