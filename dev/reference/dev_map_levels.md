# Map level of every cell of a deviation array

The map says which cells are estimated and which are fixed, and an array
with no map entry has every cell estimated. Reading it this way gives
one vector either way, `NA` at a fixed cell, which is what the dsem
checks read.

## Usage

``` r
dev_map_levels(input_list, par_name)
```

## Arguments

- input_list:

  List with `par` and `map`.

- par_name:

  Name of the deviation array.

## Value

Integer vector, one level per cell of the array, `NA` where fixed.
