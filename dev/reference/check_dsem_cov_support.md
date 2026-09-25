# Refuse covariate values outside their family's support

Refuse covariate values outside their family's support

## Usage

``` r
check_dsem_cov_support(y, family, name)
```

## Arguments

- y:

  Observations with NA for missing years.

- family:

  Family code.

- name:

  Covariate name for the message.

## Value

`invisible(NULL)`; stops on a value the family cannot produce.
