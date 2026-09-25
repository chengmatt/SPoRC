# Covariate family and link codes

The dsem package's codes and names. Families: fixed 0, gaussian 1
(normal is the same), bernoulli 2 (binomial is the same), poisson 3,
Gamma 4 (gamma is the same), gaussian_fixed_sd 5, lognormal 6, tweedie
7. Links: identity 0, log 1, logit 2, cloglog 3. `dsem_default_link`
gives each family the link its dsem constructor defaults to, except the
Gamma, whose stats default (inverse) dsem has no code for, so it gets
the log.

## Usage

``` r
dsem_family_codes()

dsem_link_codes()

dsem_default_link(family)
```

## Arguments

- family:

  Family code, for `dsem_default_link`.

## Value

Named integer vector, or one link code.
