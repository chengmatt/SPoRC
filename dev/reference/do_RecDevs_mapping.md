# Map annual recruitment deviation parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `ln_RecDevs`, the log-scale
annual recruitment deviations. The `ln_RecDevs` array has dimensions
`[n_pop x n_regions x n_years]`.

## Usage

``` r
do_RecDevs_mapping(input_list, RecDevs_spec, rec_dd)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$rec_region_prop_spec`, `$data$natal_region`, `$data$rec_dd`,
  and `$data$n_pop` to be set by upstream setup functions.

- RecDevs_spec:

  Character string specifying the sharing structure for `ln_RecDevs`, or
  `NULL` to estimate all deviations independently across all dimensions.
  Options when non-`NULL`:

  `"est_shared_r"`

  :   Separate year-specific deviations per population, shared across
      regions within each population. All regions of a given population
      follow the same annual deviation time series. Valid under both
      local and global density dependence.

  `"est_shared_pop_r"`

  :   A single set of year-specific deviations shared across all
      populations and regions. Each year receives one estimated
      parameter regardless of how many populations or regions are
      modelled. Required when `rec_dd = "global"` and `n_regions > 1`.

  `"fix"`

  :   All `ln_RecDevs` parameters fixed at zero (mapped to `NA`).
      Equivalent to deterministic recruitment with no interannual
      variability.

  `NULL`

  :   Estimate all deviations independently across populations, regions,
      and years. Not permitted when `rec_region_prop_spec = 1` and
      `n_pop > 1`, as non-natal regions have no recruitment.

- rec_dd:

  Recruitment density-dependence structure inherited from
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).
  `"global"` restricts valid `RecDevs_spec` choices to `"est_shared_r"`
  or `"est_shared_pop_r"` when `n_regions > 1`.

## Value

The input `input_list` with `$map$ln_RecDevs` set to a factor vector of
length `prod(dim(par$ln_RecDevs))`. Active parameters receive sequential
integer indices; non-natal region slots (when
`rec_region_prop_spec = 1`) and fixed deviations are `NA`. Starting
values in `$par$ln_RecDevs` are reset to `0` for any fixed cells.

## Details

Mapping behaviour is governed by two interacting considerations:

1.  **Sharing specification** (`RecDevs_spec`): controls whether
    deviations are shared across regions and/or populations, or
    estimated independently.

2.  **No-dispersal constraint**: when `rec_region_prop_spec = 1` and
    `n_pop > 1`, non-natal regions receive no recruitment and their
    deviations are structurally zero; these are automatically fixed to
    `NA` regardless of `RecDevs_spec`, and the remaining indices are
    re-numbered sequentially.

## See also

[`do_InitDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_InitDevs_mapping.md)
for the analogous initial age-structure deviation mapping, which shares
the same sharing options and no-dispersal constraint logic.
