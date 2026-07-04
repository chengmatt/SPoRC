# Simulate from a logistic-normal distribution

Draws a single composition vector of length \\K\\ from a logistic-normal
distribution using the additive log-ratio (ALR) parameterisation. The
last category is treated as the reference bin: the mean vector of the
underlying multivariate normal is \\\mu_k = \log(p_k / p_K)\\, \\k = 1,
\ldots, K-1\\. The draw is back-transformed via the additive softmax so
the returned vector sums to 1. Four covariance structures are supported,
controlled by `comp_like`.

## Usage

``` r
rlogistnormal(exp, pars, comp_like, n_sexes)
```

## Arguments

- exp:

  Numeric vector of length \\K\\. Expected (predicted) composition
  proportions; must be positive and sum to 1. The last element is used
  as the ALR reference bin and is not directly simulated.

- pars:

  Numeric vector of parameters governing the covariance structure.
  Required elements depend on `comp_like`:

  `comp_like = 2` (iid)

  :   1 parameter: `pars[1]` = marginal standard deviation \\\sigma\\.

  `comp_like = 3` (AR1 by bin)

  :   2 parameters: `pars[1]` = \\\sigma\\, `pars[2]` = AR1 correlation
      \\\rho\\ across bins. The covariance matrix is \\\sigma^2 / (1 -
      \rho^2)\\ times the AR1 correlation matrix, with the last
      row/column removed.

  `comp_like = 4` (AR1 × constant sex)

  :   3 parameters: `pars[1]` = \\\sigma\\, `pars[2]` = AR1 correlation
      across bins, `pars[3]` = constant correlation across sexes. The
      covariance is a Kronecker product \\C\_{\text{sex}} \otimes
      C\_{\text{AR1}}\\ scaled by \\\sigma^2 / ((1 -
      \rho\_{\text{bin}}^2)(1 - \rho\_{\text{sex}}^2))\\, with the last
      row/column removed.

- comp_like:

  Integer. Covariance structure: `2` = iid, `3` = AR1 by bin, `4` = AR1
  by bin with constant sex correlation.

- n_sexes:

  Integer. Number of sexes. Used only when `comp_like = 4` to construct
  the Kronecker product covariance.

## Value

Numeric vector of length \\K\\ summing to 1, representing a single draw
from the logistic-normal distribution with the specified mean and
covariance structure.
