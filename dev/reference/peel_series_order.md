# Order series so that everything needed within a year comes first

Order series so that everything needed within a year comes first

## Usage

``` r
peel_series_order(n_series, edge_from, edge_to)
```

## Arguments

- n_series:

  Number of series.

- edge_from, edge_to:

  Series indices of every dependency inside one year: the same-year
  paths, and each moderated arrow's moderating series.

## Value

List with `order`, `NULL` when the edges form a loop, and `stuck`, the
series left in that loop.
