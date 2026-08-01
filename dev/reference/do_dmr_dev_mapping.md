# Map discard mortality deviation parameters

Constructs the `logit_dmr_devs` factor map, assigning unique estimation
indices to region-year-season-fleet cells that are fished and mapping
true closures to `NA`. A cell is fished when aggregated or any
population-specific catch is used, or when the aggregate catch
observation is missing (`NA`) rather than a recorded zero, which is the
same condition under which the objective computes a non-zero `dmr`.

## Usage

``` r
do_dmr_dev_mapping(input_list, dmr_dev_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$UseCatch`, `$data$UseCatch_pop`, and `$data$ObsCatch` to be
  populated by
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- dmr_dev_spec:

  Character string specifying whether to estimate or fix deviations.
  Currently supports `"est_all"` and `"fix"`.

## Value

The input `input_list` with `$map$logit_dmr_devs` set to a factor
vector. Fished cells are assigned unique integer indices; true closures
are `NA`.

## Details

Discard observations are not required for a deviation to be estimable.
`dmr` enters the likelihood through total mortality (`ZAA`), so it is
informed by retained catch, indices, and compositions in any fished cell
where retention is less than one.
[`get_dmr_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_dmr_penalty.md)
keys on the same condition, so the estimated and penalized sets
coincide.
