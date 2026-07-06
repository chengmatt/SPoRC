# Binomial CDF, P(X \<= x), via the regularized incomplete beta

\\P(X \le x) = 1 - I_p(x + 1, n - x)\\, matching Trijoulet et al. (2023)
`dists::pbinom`. AD-safe through `RTMB::pbeta`.

## Usage

``` r
osa_pbinom(x, n, prob)
```

## Arguments

- x:

  Count.

- n:

  Number of trials.

- prob:

  Success probability.

## Value

Lower-tail CDF value.
