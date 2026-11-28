# Check that every row of a movement array sums to one

A row is where the fish in one region go, so its fractions sum to one.
Rows that do not would add or remove fish every time the matrix is
applied.

## Usage

``` r
check_movement_rows(x, arg_name)
```

## Arguments

- x:

  Movement array whose second dim is the region fish leave and whose
  third is the region they arrive in.

- arg_name:

  Character. The argument `x` was given as, for the error message.

## Value

No return value. Stops when a row does not sum to one (1e-6).
