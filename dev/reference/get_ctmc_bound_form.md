# Names for the CTMC generator bound forms

Names for the CTMC generator bound forms

## Usage

``` r
get_ctmc_bound_form(x)
```

## Arguments

- x:

  The `ctmc_diffusion_bounds` value: the code `0`, `1` or `2`, or the
  matching name `"none"`, `"softplus"` or `"upwind"`.

## Value

Character scalar naming the form, or `NA_character_` if `x` is not one
of the accepted values.
