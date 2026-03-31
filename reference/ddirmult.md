# Dirichlet Multinomial Likelihood

From
https://github.com/James-Thorson/CCSRA/blob/main/inst/executables/CCSRA_v9.cpp

## Usage

``` r
ddirmult(obs, pred, Ntotal, ln_theta, give_log = TRUE)
```

## Arguments

- obs:

  Vector of observed values in proportions

- pred:

  Vector or predicted values in proportions

- Ntotal:

  Input sample size scalar

- ln_theta:

  Weighting parameter in log space

- give_log:

  Whether or not likelihood is in log space

## Value

returns likelihood values from a dirihclet multinomial
