# Map survey length composition correlation parameters

Internal helper called by
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)
to construct the TMB/RTMB factor maps for `SrvLen_corr_pars`
`[n_regions × n_sexes × n_srv_fleets × 2]` and `SrvLen_corr_pars_agg`
`[n_srv_fleets]`, the correlation parameters for logistic-normal survey
length composition likelihoods. Follows identical activation logic to
[`do_SrvAge_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_corr_pars_mapping.md)
but operates on `SrvLenComps_LikeType`, `SrvLenComps_Type`, and
`UseSrvLenComps`.

## Usage

``` r
do_SrvLen_corr_pars_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$n_sexes`,
  `$data$SrvLenComps_Type`, `$data$SrvLenComps_LikeType`, and
  `$data$UseSrvLenComps`.

## Value

The input `input_list` with `$map$SrvLen_corr_pars` and
`$map$SrvLen_corr_pars_agg` set to factor vectors. Active parameters
receive sequential integer indices; inactive parameters are `NA`.

## See also

[`do_SrvAge_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_corr_pars_mapping.md)
for the analogous age composition correlation mapping;
[`do_SrvLen_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvLen_theta_mapping.md)
for the overdispersion parameter mapping.
