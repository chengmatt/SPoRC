# Map tag recapture overdispersion parameter

Internal helper called by
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
to construct the TMB/RTMB factor map for `ln_conv_fish_tag_theta`, the
log-scale overdispersion parameter for the tag recapture likelihood.
This is a scalar parameter active only for negative-binomial and
Dirichlet-multinomial likelihoods.

## Usage

``` r
do_conv_tag_theta_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$use_conv_fish_tagging` and `$data$conv_fish_tag_like`.

## Value

The input `input_list` with `$map$ln_conv_fish_tag_theta` set to a
length-1 factor. Active: factor level `1`; fixed: `NA`.

## Details

When no fishery fleet uses conventional tagging, the parameter is mapped
to `NA`. When tagging is active, activation depends on
`conv_fish_tag_like`: `NA` for Poisson (`0`) and multinomial (`2`, `3`);
factor level `1` for negative binomial (`1`) and Dirichlet-multinomial
(`4`, `5`).
