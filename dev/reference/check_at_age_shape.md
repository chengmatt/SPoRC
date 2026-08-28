# Refuse an at-age array or parameter that is missing a dimension

The at-age streams carry a sex margin, and nothing promotes an array
into it: an array one dimension short would otherwise be indexed by
position and read the wrong age or sex. This reports what was supplied
against what is wanted.

## Usage

``` r
check_at_age_shape(x, want, what)
```

## Arguments

- x:

  The array to check, or `NULL` to skip.

- want:

  Integer vector of the dimensions expected.

- what:

  Name used in the message.

## Value

`invisible(NULL)`. Called for its error.
