# Handle the optional modules when they were not set up

Movement, tagging and the dsem are all optional. A model with one region
has no movement to estimate and a model with no tagging data needs none
of the tagging settings, so
[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
calls those two setups with their off settings rather than requiring the
user to. Movement is only defaulted for one region, since more than one
has no population dynamics without it. The dsem needs no settings when
it was not set up, since the objective skips its density entirely, so
the only case to catch is a module that declared `"dsem"` without one
being given.

## Usage

``` r
fill_optional_modules(data, parameters, mapping)
```

## Arguments

- data, parameters, mapping:

  The three lists
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
  was given.

## Value

List with `data`, `par` and `map`, unchanged where the module was
already set up.
