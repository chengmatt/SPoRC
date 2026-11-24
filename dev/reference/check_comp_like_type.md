# Check a composition likelihood setting

Names the fleets that were given something unrecognized, and suggests
the accepted value when the entry differs only in case, spacing or
punctuation.

## Usage

``` r
check_comp_like_type(
  x,
  what,
  allowed = c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal",
    "1d-Logistic-Normal", "2d-Logistic-Normal", "iid-Logistic-Normal-miss0",
    "1d-Logistic-Normal-miss0", "2d-Logistic-Normal-miss0"),
  note =
    paste0("The -miss0 forms drop the empty bins and scale the standard deviation by the ",
    "input sample size.")
)
```

## Arguments

- x:

  Character vector given for the setting, one entry per fleet.

- what:

  Character. The argument name, used in the message.

- allowed:

  Character vector of accepted values, in the order the message lists
  them.

- note:

  Character. Sentence appended after the accepted values, explaining
  what the less obvious ones do.

## Value

`x`, invisibly.
