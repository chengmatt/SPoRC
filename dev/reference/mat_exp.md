# Matrix exponential of a movement generator, exactly or by implicit solve

Turns an instantaneous rate matrix into transition fractions. This is
the single point where `SPoRC` decides how a matrix exponential is
evaluated, so that
[`Get_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Movement.md)
and every `move_timing = 2` operator in `model_transition.R` share one
convention.

## Usage

``` r
mat_exp(A, expm_nsub = 0)
```

## Arguments

- A:

  Square matrix, dense or sparse, numeric or `advector`. The generator
  whose exponential is wanted, already scaled by whatever time step the
  caller intends (i.e. this returns \\e^{A}\\, not \\e^{A\Delta}\\).

- expm_nsub:

  Integer. `0` (default) evaluates \\e^{A}\\ with
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html). A
  power of two \\n \ge 1\\ uses the implicit (backward Euler) scheme
  with \\n\\ substeps (faster but loses acurracy).

## Value

A plain dense matrix of the same dimension as `A`.
