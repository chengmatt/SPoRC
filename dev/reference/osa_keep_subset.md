# Keep-subset for internal OSA residuals

Elements flagged `last_in_group == TRUE` are the statistically
determined/reference cell of their group (excluded from the discrete OSA
evaluation).

## Usage

``` r
osa_keep_subset(last_in_group)
```

## Arguments

- last_in_group:

  Logical vector (with possible NAs), as produced by
  `pack_comp_osa(..., return_labels = TRUE)` or
  `pack_tag_osa(..., return_labels = TRUE)`.

## Value

Integer vector of positions to keep.
