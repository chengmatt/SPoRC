# Map recruitment regional apportionment parameters

Builds the factor map for `rec_region_prop_pars`
`[n_pop x (n_regions - 1)]`, the logit-scale share of recruits per
region, under a softmax with one reference region omitted. Both the
parameter and its map are `NULL` when `n_regions = 1`. Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

## Usage

``` r
do_rec_region_prop_mapping(input_list, rec_region_prop_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`. Requires `$data$n_pop`,
  `$data$n_regions` and `$data$natal_region`.

- rec_region_prop_spec:

  Dispersal structure. `"no_dispersal"` assigns recruits entirely to
  their natal region, overwriting the starting values with `-20` for
  non-natal regions and `+20` for the natal region when
  `natal_region > 1`, and mapping every element to `NA`. `NULL`
  estimates all independently. Both require `n_regions > 1`.

## Value

`input_list` with `$map$rec_region_prop_pars` set to a factor vector of
length `n_pop * (n_regions - 1)`, or `NULL` when `n_regions = 1`.
Starting values are overwritten under `"no_dispersal"`.
