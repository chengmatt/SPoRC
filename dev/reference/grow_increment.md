# Grow a mean length forward by a fraction of a year

The increment of the Richards (or von Bertalanffy, `rho = 1`) curve over
an elapsed time `e` from a current mean length `L`: \$\$L(t + e)^\rho =
L\_\infty^\rho + (L^\rho - L\_\infty^\rho) e^{-K e}\$\$ Applied to a
length that sits on the curve it returns the curve's value `e` later, so
splitting a year into seasons changes nothing; applied to a length kept
from an earlier year's parameters it is how a cohort keeps its own
history.

## Usage

``` r
grow_increment(L, e, K, Linf, rho = 1)
```

## Arguments

- L:

  Mean length(s) at the start, possibly AD.

- e:

  Elapsed time in years (data).

- K, Linf, rho:

  Growth rate, asymptote and Richards coefficient in effect.

## Value

Mean length(s) after the increment.
