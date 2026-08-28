# Penalty on the level of the recruitment series itself

Penalizes the log recruitment time series directly, independently of any
stock-recruit residual penalty. Under a stock-recruit relationship the
deviations are residuals about the predicted curve, so a model that also
wants to keep the recruitment series itself from wandering has nowhere
to say so; this is that second statement. Centering on the series' own
mean penalizes only its variability and leaves its level to the rest of
the model.

## Usage

``` r
get_rec_level_penalty(Rec, sigma, center = 1, yrs = NULL)
```

## Arguments

- Rec:

  Array `[pop, region, year]` of recruitment.

- sigma:

  Numeric standard deviation of the penalty.

- center:

  Integer. `1` centers on the mean of the log series, `0` centers on
  zero.

- yrs:

  Integer vector of years the penalty applies over, or `NULL` for every
  year.

## Value

Array `[pop, region, year]` of negative log-likelihood contributions,
zero outside `yrs`.
