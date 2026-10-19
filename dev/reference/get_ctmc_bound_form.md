# Names for the CTMC generator bound forms

Names for the CTMC generator bound forms

## Usage

``` r
get_ctmc_bound_form(x)
```

## Arguments

- x:

  The `ctmc_diffusion_bounds` value: the code `0`, `1`, `2` or `3`, or
  the matching name `"none"`, `"softplus"`, `"upwind"` or `"sg"`.

## Value

Character scalar naming the form, or `NA` if `x` is not one of the
accepted values.
