# Map initial tag-induced mortality parameter

Internal helper called by
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
to construct the TMB/RTMB factor map for `ln_init_conv_tag_mort`, the
log-scale mortality applied to fish at the moment of tag release. This
is a scalar parameter (not dimensioned by fleet or region).

## Usage

``` r
do_conv_init_tag_mort_mapping(input_list, init_conv_tag_mort_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$use_conv_fish_tagging`.

- init_conv_tag_mort_spec:

  Character string. One of:

  `"fix"`

  :   Fix `ln_init_conv_tag_mort` at its starting value (mapped to
      `NA`).

  `"est"`

  :   Estimate `ln_init_conv_tag_mort` (mapped to factor level `1`).

## Value

The input `input_list` with `$map$ln_init_conv_tag_mort` set to a
length-1 factor. Active: factor level `1`; fixed: `NA`.

## Details

When no fishery fleet uses conventional tagging
(`all(use_conv_fish_tagging == 0)`), the parameter is automatically
mapped to `NA` regardless of `init_conv_tag_mort_spec`.
