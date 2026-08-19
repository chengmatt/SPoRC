# Stock-recruit residual penalty under mean recruitment

Compares the recruitment series against a stock-recruit curve without
letting the curve generate it. Under `rec_model = "mean_rec"`
recruitment is \\R_y = \exp(\mu + \varepsilon_y)\\ and the curve enters
only here, as a Gaussian on the log residual \\\log R_y - \log
\widehat{R}\_y\\. That is a different statement from the deviation
penalty applied when the curve generates recruitment: there the residual
is the parameter, here it is a derived quantity and the deviations
remain free.

## Usage

``` r
get_sr_penalty(Rec, SR_pred, sigma, yrs = NULL)
```

## Arguments

- Rec:

  Array `[pop, region, year]` of realized recruitment.

- SR_pred:

  Array `[pop, region, year]` of the curve's prediction, computed
  alongside the population projection.

- sigma:

  Numeric standard deviation of the residual.

- yrs:

  Integer vector of years the penalty applies over, or `NULL` for every
  year. Years outside it contribute zero and stay free.

## Value

Array `[pop, region, year]` of negative log-likelihood contributions,
zero outside `yrs`.

## Details

Several AFSC models are written this way to reflect that a weakly
determined SR relationship
