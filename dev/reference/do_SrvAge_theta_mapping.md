# Map survey age composition overdispersion parameters

Internal helper called by
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)
to construct the TMB/RTMB factor maps for `ln_SrvAge_theta`
`[n_regions × n_sexes × n_srv_fleets]` and `ln_SrvAge_theta_agg`
`[n_srv_fleets]`, the log-scale overdispersion parameters for survey age
composition likelihoods.

## Usage

``` r
do_SrvAge_theta_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$n_sexes`,
  `$data$SrvAgeComps_Type`, `$data$SrvAgeComps_LikeType`, and
  `$data$UseSrvAgeComps`.

## Value

The input `input_list` with `$map$ln_SrvAge_theta` and
`$map$ln_SrvAge_theta_agg` set to factor vectors. Active parameters
receive sequential integer indices; inactive parameters are `NA`.

## Details

Parameters are activated only when the fleet's likelihood type is not
multinomial (`SrvAgeComps_LikeType != 0`) and at least one year of age
composition data is used (`UseSrvAgeComps > 0`). Within active fleets,
the composition type (`SrvAgeComps_Type`) determines which parameter
array is populated: aggregated compositions use `ln_SrvAge_theta_agg`;
split or joint-by-sex compositions use `ln_SrvAge_theta`. Fleets using
multinomial likelihoods or with no active age composition data have all
parameters mapped to `NA`.

## See also

[`do_SrvLen_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvLen_theta_mapping.md)
for the analogous length composition overdispersion mapping;
[`do_SrvAge_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_corr_pars_mapping.md)
for the associated correlation parameter mapping.
