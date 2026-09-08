# Share of one penalty owed by each cell of a mapped deviation array

Cells sharing a map level hold one parameter between them, so penalizing
every cell would count that parameter once per cell. Each cell takes one
over the number of cells at its level, and a cell mapped off takes zero.

## Usage

``` r
dev_share_weights(map, dims)
```

## Arguments

- map:

  Numeric array of map levels, `NA` where the cell is fixed, or `NULL`
  when every cell holds its own parameter.

- dims:

  Dimensions the weights are returned on.

## Value

Numeric array of weights over `dims`, all ones when `map` is `NULL`.
