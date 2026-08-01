# Convert character or numeric input to numeric codes

Maps a character vector to integer codes via a named lookup, or passes
numeric input through unchanged. Character arrays and matrices are
converted element-wise and keep their original dimensions. Unrecognised
character values raise an informative error listing both the invalid
inputs and the valid options.

## Usage

``` r
convert_to_numeric(x, lookup)
```

## Arguments

- x:

  Character vector, numeric vector, or array to convert.

- lookup:

  Named list (or named atomic vector) mapping valid character strings to
  numeric codes (e.g., `list("none" = 999, "multinomial" = 0)`).

## Value

Numeric vector or array of the same shape as `x`, without names.
