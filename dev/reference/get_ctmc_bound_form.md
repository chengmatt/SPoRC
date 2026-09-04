# Names for the CTMC generator bound forms

`ctmc_diffusion_bounds` started life as a 0/1 flag and now also takes
the name of the form directly, so input lists built before the named
forms existed keep evaluating. Shared by `Setup_Mod_Movement` and
`Get_Movement` so validation and the objective agree on the spelling.

## Usage

``` r
get_ctmc_bound_form(x)
```

## Arguments

- x:

  The `ctmc_diffusion_bounds` value: `0`, `1`, or one of `"none"`,
  `"softplus"`, `"clamp"`, `"upwind"`, `"barker"`, `"logsoftplus"`.

## Value

Character scalar naming the form, or `NA_character_` if `x` is not one
of the accepted values.
