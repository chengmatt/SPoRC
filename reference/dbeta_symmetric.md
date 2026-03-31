# Symmetric Beta Function

Symmetric Beta Function

## Usage

``` r
dbeta_symmetric(p_val, p_ub, p_lb, p_prsd, log = TRUE)
```

## Arguments

- p_val:

  Parameter value

- p_ub:

  Upper Bound of Parameter

- p_lb:

  Lower Bound of Parameter

- p_prsd:

  SD of parameter, higher values have a stronger penalty on bounds,
  lower values have a more difuse penalty on bounds

- log:

  whether or not to return the log likelihood

## Value

Returns likelihood values from a symmetric beta distribution
