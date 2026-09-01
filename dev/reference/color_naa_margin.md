# Apply a correlation factor along one margin of an array

Colors an array of independent normals so that one margin carries a
given correlation, the reverse of the whitening the penalty uses.
Operating on plain doubles rather than on the AD tape, so the margin can
simply be permuted to the front rather than being reached by index
arithmetic.

## Usage

``` r
color_naa_margin(x, L, margin)
```

## Arguments

- x:

  Numeric array.

- L:

  Lower triangular factor of the margin's correlation matrix.

- margin:

  Integer margin to apply it along.

## Value

An array of the same shape.
