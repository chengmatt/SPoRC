# The sd line of each linked series

A series' own sd is the two-headed arrow from it to itself at lag zero.
Under `RecDevs_model = "dsem"` that value stands in for `sigmaR`, so the
objective needs to know which arrow it is. A moderated sd changes by
year and cannot stand in, so that series is given a 0.

## Usage

``` r
get_dsem_link_sd_arrow(dsem_model, link_col)
```

## Arguments

- dsem_model:

  From
  [`read_dsem_arrows`](https://chengmatt.github.io/SPoRC/dev/reference/read_dsem_arrows.md).

- link_col:

  Grid column of each linked series.

## Value

Integer vector, one arrow row per linked series, 0 where none serves.
