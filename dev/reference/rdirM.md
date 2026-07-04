# Draw samples from a Dirichlet-multinomial distribution

Generates `n` independent draws from a Dirichlet-multinomial
distribution with total count `N` and concentration parameter vector
`alpha`. Each draw is produced by first sampling a Dirichlet-
distributed probability vector from `Gamma(alpha)` variates and then
drawing multinomial counts conditioned on that probability vector.

## Usage

``` r
rdirM(n, N, alpha)
```

## Arguments

- n:

  Integer. Number of independent draws to generate.

- N:

  Integer. Total count per draw (multinomial sample size).

- alpha:

  Numeric vector of length \\K\\. Dirichlet concentration parameters.
  All elements must be positive. Larger values relative to \\N\\ produce
  draws closer to the expected proportions \\\alpha / \sum \alpha\\.

## Value

Integer matrix of dimensions \\K \times n\\. Each column is one draw: a
vector of category counts summing to `N`.

## Details

The Dirichlet-multinomial arises naturally as an overdispersed
alternative to the multinomial: the marginal variance of each category
count is \\N \bar{p}\_k (1 - \bar{p}\_k) (N + \theta) / (1 + \theta)\\,
where \\\theta = \sum \alpha_k\\ controls overdispersion. In SPoRC,
`alpha` is typically supplied as \\\exp(\ln\theta) \times N \times
\hat{p}\\ so that \\\exp(\ln\theta)\\ is the per-observation
overdispersion scalar.
