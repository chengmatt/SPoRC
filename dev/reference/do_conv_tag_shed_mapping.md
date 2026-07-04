# Map chronic tag shedding parameter

Internal helper called by
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
to construct the TMB/RTMB factor map for `ln_conv_tag_shed`, the
log-scale annual chronic tag shedding rate. This is a scalar parameter.

## Usage

``` r
do_conv_tag_shed_mapping(input_list, conv_tag_shed_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$use_conv_fish_tagging`.

- conv_tag_shed_spec:

  Character string. One of:

  `"fix"`

  :   Fix `ln_conv_tag_shed` at its starting value (mapped to `NA`).

  `"est"`

  :   Estimate `ln_conv_tag_shed` (mapped to factor level `1`).

## Value

The input `input_list` with `$map$ln_conv_tag_shed` set to a length-1
factor. Active: factor level `1`; fixed: `NA`.

## Details

When any fishery fleet has tagging disabled
(`any(use_conv_fish_tagging == 0)`), the parameter is automatically
mapped to `NA` regardless of `conv_tag_shed_spec`.
