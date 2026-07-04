# Map survey length composition overdispersion parameters

Internal helper called by
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)
to construct the TMB/RTMB factor maps for `ln_SrvLen_theta`
`[n_regions × n_sexes × n_srv_fleets]` and `ln_SrvLen_theta_agg`
`[n_srv_fleets]`, the log-scale overdispersion parameters for survey
length composition likelihoods. Follows identical activation logic to
[`do_SrvAge_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_theta_mapping.md)
but operates on `SrvLenComps_LikeType`, `SrvLenComps_Type`, and
`UseSrvLenComps`.

## Usage

``` r
do_SrvLen_theta_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$n_sexes`,
  `$data$SrvLenComps_Type`, `$data$SrvLenComps_LikeType`, and
  `$data$UseSrvLenComps`.

## Value

The input `input_list` with `$map$ln_SrvLen_theta` and
`$map$ln_SrvLen_theta_agg` set to factor vectors. Active parameters
receive sequential integer indices; inactive parameters are `NA`.

## See also

[`do_SrvAge_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_theta_mapping.md)
for the analogous age composition overdispersion mapping;
[`do_SrvLen_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvLen_corr_pars_mapping.md)
for the associated correlation parameter mapping.
