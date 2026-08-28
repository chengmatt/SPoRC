# Do a fleet's at-age observations fill a complete age by year grid?

A separable correlation over ages and years is a statement about a
rectangle. This reports whether every combination of the years and ages
a fleet observes is present, within each population, region, season and
sex.

## Usage

``` r
at_age_block_complete(use_arr, f, nd, i_y, i_a)
```

## Arguments

- use_arr:

  Use array for the stream.

- f:

  Fleet index.

- nd:

  Number of dimensions of `use_arr`.

- i_y, i_a:

  Positions of the year and age dimensions.

## Value

`TRUE` when every block is complete.
