# Draw the dsem grid year by year

A moderated arrow puts the grid's own values inside \\B\\, so there is
no one covariance to draw from. Each year is drawn instead in the order
`series_order` gives, every cell from what points into it and an
innovation, which is the model the density evaluates. Years already
known are held as is.

## Usage

``` r
draw_dsem_recursive(
  dsem_model,
  arrow_value,
  mu_grid,
  n_sims,
  n_cond = 0,
  x_known = NULL
)
```

## Arguments

- dsem_model:

  Output of `read_dsem_arrows`.

- arrow_value:

  Values of the arrows, from `get_dsem_arrow_values`.

- mu_grid:

  Matrix `[year, series]` of series means.

- n_sims:

  Number of replicates.

- n_cond:

  Leading years kept at `x_known` rather than drawn.

- x_known:

  Matrix `[year, series]` read for those years.

## Value

Array `[year, series, sim]`.
