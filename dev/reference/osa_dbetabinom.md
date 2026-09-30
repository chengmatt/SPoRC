# Two-category Dirichlet-multinomial (beta-binomial) log-density

Taken on the number of trials rather than on a pair of counts, so that
the trials stay fixed while
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
sweeps `x` over its support. A candidate above `size` carries no mass,
and `-lgamma(rest + 1)` is what says so: it is the one term with a pole
there, so the shape term is held at zero rather than allowed to answer
with a second pole that would cancel it and leave `NaN` behind.

## Usage

``` r
osa_dbetabinom(x, size, shape1, shape2)
```

## Arguments

- x:

  Count in the first category.

- size:

  Number of trials.

- shape1, shape2:

  Concentration of the first category and of the rest.

## Value

Scalar log-density.
