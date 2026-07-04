# Convert character or numeric input to numeric codes

Maps a character vector to integer codes via a named lookup list, or
passes numeric input through unchanged. Arrays and matrices are
flattened, converted element-wise, and restored to their original
dimensions. Unrecognised character values raise an informative error
listing both the invalid inputs and the valid options.

## Usage

``` r
convert_to_numeric(x, lookup)
```

## Arguments

- x:

  Character vector, numeric vector, or array to convert.

- lookup:

  Named list mapping valid character strings to numeric codes (e.g.,
  `list("none" = 999, "multinomial" = 0)`).

## Value

Numeric vector or array of the same shape as `x`.
