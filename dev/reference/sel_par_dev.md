# A selectivity parameter in one year

The parameter times `exp(dev)`, or the parameter plus `dev` where the
fleet asks for deviations added on the natural scale.

## Usage

``` r
sel_par_dev(value, dev, additive)
```

## Arguments

- value:

  Parameter on its natural scale, for example an age at 50% selection.

- dev:

  That year's deviation.

- additive:

  `TRUE` adds the deviation, `FALSE` multiplies by its exponential.

## Value

The parameter's value that year.
