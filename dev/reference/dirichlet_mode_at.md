# A Dirichlet moved so its mode sits at given fractions

The concentration (the sum of `alpha`) is kept, so the prior keeps its
strength and loses its pull at `frac`. A Dirichlet whose concentration
is at or below its number of cells has no interior mode and is returned
as given.

## Usage

``` r
dirichlet_mode_at(alpha, frac)
```

## Arguments

- alpha:

  Concentration vector.

- frac:

  Fractions summing to one, the same length.

## Value

The moved `alpha`.
