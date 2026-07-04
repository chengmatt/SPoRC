# Two-category Dirichlet-multinomial log-density

Computes the log-probability mass function of a two-category
Dirichlet-multinomial (beta-binomial) distribution. AD-safe and intended
for use as the conditional building block of
[`ddirmult_osa`](https://chengmatt.github.io/SPoRC/dev/reference/ddirmult_osa.md).

## Usage

``` r
ddirmult2(obs2, alpha2)
```

## Arguments

- obs2:

  Numeric/AD vector of length 2 with observed counts
  `c(count_a, count_remaining)`.

- alpha2:

  Numeric/AD vector of length 2 with Dirichlet concentration parameters
  `c(alpha_a, alpha_remaining)`.

## Value

Scalar log-likelihood contribution.

## Details

The closed-form log-pmf is \$\$\log\Gamma(N + 1) - \sum_i
\log\Gamma(x_i + 1) + \log\Gamma(A_0) - \log\Gamma(N + A_0) + \sum_i
\left\[ \log\Gamma(x_i + \alpha_i) - \log\Gamma(\alpha_i) \right\],\$\$
with \\N = \sum_i x_i\\ and \\A_0 = \sum_i \alpha_i\\.
