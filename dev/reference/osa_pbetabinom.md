# Beta-binomial CDF, P(X \<= x), by summation

Sums two-category Dirichlet-multinomial (beta-binomial) probabilities
over \\0, \ldots, x\\, matching Trijoulet et al. (2023) `pbetabinom`.
Clamped at 1 to guard against floating-point overshoot.

## Usage

``` r
osa_pbetabinom(x, N, alpha, beta)
```

## Arguments

- x:

  Count (upper summation limit; taken from frozen observed values).

- N:

  Number of trials.

- alpha:

  First beta-binomial shape.

- beta:

  Second beta-binomial shape.

## Value

Lower-tail CDF value.
