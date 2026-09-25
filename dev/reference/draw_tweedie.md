# Draw from a tweedie as a poisson sum of gammas

A tweedie with power between 1 and 2 is a poisson number of gamma jumps:
the count has rate \\\mu^{2-p} / (\phi (2-p))\\, each jump shape
\\(2-p)/(p-1)\\ and scale \\\phi (p-1) \mu^{p-1}\\, and no jumps is a
zero. Mean \\\mu\\, variance \\\phi \mu^p\\.

## Usage

``` r
draw_tweedie(mu, phi, p)
```

## Arguments

- mu:

  Means, one per draw.

- phi:

  Dispersion.

- p:

  Power in (1, 2).

## Value

Numeric vector, one draw per mean.
