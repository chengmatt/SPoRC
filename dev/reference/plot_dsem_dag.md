# Plot a dsem as a graph of its arrows

One node per series, and a further node `lag(series, k)` for a source
read `k` years back. A path is a solid edge, blue when its estimate is
positive and red when negative; an sd line is a grey loop on its own
node and a covariance a grey dashed line with no head. Nodes are laid
out in layers, in base graphics.

## Usage

``` r
plot_dsem_dag(
  x,
  sd_rep = NULL,
  edge_label = c("value_and_stars", "value", "name"),
  digits = 2,
  ...
)
```

## Arguments

- x:

  Output of `read_dsem_arrows`, or a fitted model whose data holds one
  (`Setup_Mod_DSEM` then `fit_model`).

- sd_rep:

  Optional `RTMB::sdreport` of that fit, for the p values. Read from
  `x$sd_rep` when the fit has one.

- edge_label:

  What to write on each edge: the estimate with stars for its p value
  (`***` below 0.001, `**` below 0.01, `*` below 0.05), the estimate
  alone, or the parameter name. A fixed arrow shows its value in every
  case, and a bare arrow set has no estimates, so it is labeled with
  names.

- digits:

  Digits to show for an estimate.

- ...:

  Passed to
  [`igraph::plot.igraph`](https://r.igraph.org/reference/plot.igraph.html),
  such as `vertex.size`.

## Value

The `igraph` object, invisibly. Its edge attributes hold each arrow's
lag, estimate, standard error and p value.

## See also

[`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md)
