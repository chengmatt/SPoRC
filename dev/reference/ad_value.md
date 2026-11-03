# Turn Model Quantities Into Plain Numbers

Inside the model, RTMB records each calculation so it can differentiate
it, and a recorded (AD) quantity cannot be compared with `if`. A harvest
rule or a constraint that switches on stock status or catch reads a
plain copy through this function and compares that instead. The copy has
no derivative, so the outcome of the comparison is fixed, while the
quantity itself is unchanged and still differentiated wherever the rule
uses it in arithmetic. RTMB records the model once, when it is built, so
the record keeps whichever side of the comparison was true then. With a
projection set up by
[`Setup_Mod_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md),
[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
builds the model again at the estimates, so that side is the estimated
stock's; in a model of your own, build it again at the estimates
yourself. Lists are read element by element, and plain numbers are
returned unchanged.

## Usage

``` r
ad_value(x)
```

## Arguments

- x:

  A model quantity (a number, vector or array), a list of them, or any
  plain R object.

## Value

The same shape, names and dimnames, with plain numbers in place of model
quantities.

## Examples

``` r
ad_value(list(a = 1:3, b = NULL)) # unchanged outside the model
#> $a
#> [1] 1 2 3
#> 
#> $b
#> NULL
#> 
```
