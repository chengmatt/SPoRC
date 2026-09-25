# Values of the arrows from the parameter vectors

Each arrow line becomes one number here. A path or covariance line reads
its coefficient from `dsem_beta`, an sd line reads `exp(ln_dsem_sd)`, a
line written with `NA` for a name keeps the value it was given, and a
line named after a series is left at zero because its value changes year
to year and is read off the grid where the matrices are filled.

## Usage

``` r
get_dsem_arrow_values(dsem_beta, ln_dsem_sd, dsem_model)
```

## Arguments

- dsem_beta:

  Numeric vector, one entry per estimated path or covariance.

- ln_dsem_sd:

  Numeric vector, one entry per estimated sd line, log scale.

- dsem_model:

  Output of
  [`read_dsem_arrows`](https://chengmatt.github.io/SPoRC/dev/reference/read_dsem_arrows.md).

## Value

Numeric vector with one value per arrow, in the arrows' order.

## Details

The lines `env -> rec, 0, b` (`dsem_beta[1] = 0.4`), `env <-> env, 0, s`
(`ln_dsem_sd[1] = log(0.9)`) and `rec <-> rec, 0, NA, 1` give the values
0.4, 0.9 and 1.
