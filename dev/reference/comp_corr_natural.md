# Put a composition correlation parameter on its natural scale

The logistic normal forms that keep every bin take a correlation
anywhere in \\(-1, 1)\\; the forms that drop the zeros take a positive
one, matching the autoregression they mirror. One place decides which,
so an operating model draws on the same scale the estimation model fits
on.

## Usage

``` r
comp_corr_natural(pars, comp_like)
```

## Arguments

- pars:

  Unconstrained correlation parameters.

- comp_like:

  Integer likelihood code for the fleet.

## Value

The correlations on the natural scale.
