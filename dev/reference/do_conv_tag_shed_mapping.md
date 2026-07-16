# Map chronic tag shedding parameter

Internal helper called by
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
to construct the TMB/RTMB factor map for `ln_conv_tag_shed`, the
log-scale annual chronic tag shedding rate. This parameter is a vector
of length `n_conv_tag_cohorts` (one element per tag release event / row
of `conv_tag_release_indicator`); when tagging is inactive it is a
length-1 placeholder.

## Usage

``` r
do_conv_tag_shed_mapping(input_list, conv_tag_shed_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$use_conv_fish_tagging` and `$par$ln_conv_tag_shed`.

- conv_tag_shed_spec:

  Character string. One of:

  `"fix"`

  :   Fix `ln_conv_tag_shed` at its starting values for every release
      event (mapped to `NA`).

  `"est_shared"`

  :   Estimate a single value of `ln_conv_tag_shed` shared across all
      release events (mapped to factor level `1` for every event).

  `"est_all"`

  :   Estimate an independent value of `ln_conv_tag_shed` for every
      release event (mapped to distinct factor levels
      `1:n_conv_tag_cohorts`). Issues a warning, as per-event chronic
      shedding is frequently non-identifiable.

## Value

The input `input_list` with `$map$ln_conv_tag_shed` set to a factor
vector the same length as `$par$ln_conv_tag_shed`. Fixed: all `NA`;
shared: all factor level `1`; independent: distinct levels per event.

## Details

When no fishery fleet uses conventional tagging
(`all(use_conv_fish_tagging == 0)`), the parameter is automatically
mapped to `NA` regardless of `conv_tag_shed_spec`.
