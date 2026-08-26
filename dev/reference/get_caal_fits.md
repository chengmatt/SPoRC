# Conditional age-at-length fits as mean age within each length bin

A conditional age-at-length observation is an age composition within one
length bin, and a model can carry thousands of them, so the composition
itself is not something to plot row by row. The mean age within each
length bin is the summary that reads: it is what the data say about
growth, and it is the statistic Francis reweighting penalizes the fit
against. Returned observed and predicted, so the same frame draws mean
age over years and mean age against length pooled over years. Residuals
come from
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
with `comp_source = "Fish_caal"` or `"Srv_caal"`, plotted with
[`plot_resids`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md).

## Usage

``` r
get_caal_fits(data, rep)
```

## Arguments

- data:

  Data list from the fitted model.

- rep:

  Report list from the fitted model.

## Value

A data frame with one row per region, year, season, length bin, sex and
fleet that aged fish, holding `obs` and `pred` mean age, the predicted
standard deviation of age within the bin (`sd_pred`) and the number aged
(`ISS`). Empty when the model carries no conditional age-at-length data.
