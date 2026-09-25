# Map annual recruitment deviation parameters

Builds the factor map for `ln_RecDevs` `[n_pop x n_regions x n_years]`,
the log-scale annual recruitment deviations. When
`rec_region_prop_spec = 1` and `n_pop > 1`, non-natal regions get no
recruitment, so their deviations are fixed to `NA` whatever
`RecDevs_spec` asks for and the rest are re-numbered. Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

## Usage

``` r
do_RecDevs_mapping(input_list, RecDevs_spec, rec_dd, dont_pen_recdev_first = 0)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`. Requires
  `$data$rec_region_prop_spec`, `$data$natal_region`, `$data$rec_dd` and
  `$data$n_pop`.

- RecDevs_spec:

  Sharing structure for `ln_RecDevs`. `"est_shared_r"` gives one
  deviation series per population, shared across its regions.
  `"est_shared_pop_r"` gives one series across every population and
  region, required when `rec_dd = "global"` and `n_regions > 1`. `"fix"`
  holds every deviation at zero. `NULL` estimates all independently,
  which is not permitted when `rec_region_prop_spec = 1` and
  `n_pop > 1`.

- rec_dd:

  Density dependence inherited from
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).
  `"global"` restricts `RecDevs_spec` to `"est_shared_r"` or
  `"est_shared_pop_r"` when `n_regions > 1`.

## Value

`input_list` with `$map$ln_RecDevs` set to a factor vector of length
`prod(dim(par$ln_RecDevs))`. Active parameters take sequential integers;
non-natal region slots and fixed deviations are `NA`, and their starting
values are reset to `0`.

## See also

[`do_InitDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_InitDevs_mapping.md),
which shares the same options.
