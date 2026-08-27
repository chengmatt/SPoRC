# A bin selection array, or NULL when it restricts nothing

The composition machinery treats `NULL` as "fit every bin", which lets
the likelihood and the OSA packers skip the restriction entirely and
lets a backwards-compatible all-ones array of the wrong length never be
indexed into. Both the objective and
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
decide that here, so they cannot disagree about which bins were fitted.

## Usage

``` r
bins_or_null(x)
```

## Arguments

- x:

  A `[n_obs_bins x n_fleets]` 0/1 array, or `NULL`.

## Value

`x`, or `NULL` if it is absent or selects every bin.
