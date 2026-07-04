# Map recruitment regional apportionment parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `rec_region_prop_pars`, the
logit-scale parameters controlling the proportion of recruits assigned
to each region. The array has dimensions `[n_pop x (n_regions - 1)]`,
using a sum-to-one soft-max parameterisation with one reference region
omitted.

## Usage

``` r
do_rec_region_prop_mapping(input_list, rec_region_prop_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_pop`, `$data$n_regions`, and `$data$natal_region` to be set
  by upstream setup functions.

- rec_region_prop_spec:

  Character string specifying the dispersal structure, or `NULL` to
  estimate all regional proportions freely. Options when non-`NULL`:

  `"no_dispersal"`

  :   Recruits are assigned entirely to their natal region. Starting
      values for `rec_region_prop_pars` are overwritten with
      large-magnitude values (`-20` for non-natal regions, `+20` for the
      natal region when `natal_region > 1`), and all elements are mapped
      to `NA` so the parameters are not estimated. Requires `n_pop > 1`
      and `n_regions > 1`.

  `NULL`

  :   All `rec_region_prop_pars` are estimated independently. Only
      available when `n_regions > 1`.

## Value

The input `input_list` with `$map$rec_region_prop_pars` set to a factor
vector of length `n_pop * (n_regions - 1)`, or `NULL` when
`n_regions = 1`. Under `"no_dispersal"`, `$par$rec_region_prop_pars`
starting values are also overwritten.

## Details

When `n_regions = 1`, the parameter is structurally irrelevant and both
`$par$rec_region_prop_pars` and `$map$rec_region_prop_pars` are set to
`NULL`.
