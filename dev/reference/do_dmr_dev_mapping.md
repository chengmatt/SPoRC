# Map discard mortality deviation parameters

Constructs the `logit_dmr_devs` factor map, assigning unique estimation
indices to region–year–season–fleet cells where discard data are used
(`UseDiscard == 1`) and mapping cells without discard data to `NA`. This
ensures that `logit_dmr_devs` parameters are only estimated for
dimensions with observed discard.

## Usage

``` r
do_dmr_dev_mapping(input_list, dmr_dev_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$UseDiscard` and `$data$UseDiscard_pop` to be populated by
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- dmr_dev_spec:

  Character string specifying whether to estimate or fix deviations.
  Currently supports `"est_all"` and `"fix"`.

## Value

The input `input_list` with `$map$logit_dmr_devs` set to a factor
vector. Cells with discard are assigned sequential integer indices;
cells without discard are `NA`.
