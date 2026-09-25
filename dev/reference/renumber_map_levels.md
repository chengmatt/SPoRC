# Renumber a map so its levels run from one with no gaps

Blanking cells of a shared map leaves holes in the level numbering,
which `factor` keeps as unused levels.

## Usage

``` r
renumber_map_levels(x)
```

## Arguments

- x:

  Integer array holding map levels and `NA`.

## Value

The same array with its levels renumbered.
