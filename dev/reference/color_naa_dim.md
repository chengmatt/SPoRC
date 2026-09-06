# Apply a correlation factor along one dim of an array

Colors an array of independent normals so that one dim has a given
correlation, the reverse of the whitening the penalty uses. Operating on
plain doubles rather than on the AD tape, so the dim can simply be
permuted to the front rather than being reached by index arithmetic.

## Usage

``` r
color_naa_dim(x, L, dim_idx)
```

## Arguments

- x:

  Numeric array.

- L:

  Lower triangular factor of the dim's correlation matrix.

- dim_idx:

  Integer dim to apply it along.

## Value

An array of the same shape.
