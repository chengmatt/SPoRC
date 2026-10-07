# One sum of the written out multinomial

`sum(obs_w * log(props + const))`, with a bin observed empty adding zero
rather than `0 * log(0)`, which is `NaN` under `addtocomp = 0`.

## Usage

``` r
comp_mltnml_term(obs_w, props, const)
```

## Arguments

- obs_w:

  Observed proportions weighting the sum, plus `const` or not.

- props:

  Expected or observed proportions the logarithm is taken of.

- const:

  Composition constant, `addtocomp`.

## Value

Numeric scalar.
