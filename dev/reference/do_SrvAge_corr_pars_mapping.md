# Map survey age composition correlation parameters

Internal helper called by
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)
to construct the TMB/RTMB factor maps for `SrvAge_corr_pars`
`[n_regions × n_sexes × n_srv_fleets × 2]` and `SrvAge_corr_pars_agg`
`[n_srv_fleets]`, the correlation parameters for logistic-normal survey
age composition likelihoods.

## Usage

``` r
do_SrvAge_corr_pars_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$n_sexes`,
  `$data$SrvAgeComps_Type`, `$data$SrvAgeComps_LikeType`, and
  `$data$UseSrvAgeComps`.

## Value

The input `input_list` with `$map$SrvAge_corr_pars` and
`$map$SrvAge_corr_pars_agg` set to factor vectors. Active parameters
receive sequential integer indices; inactive parameters are `NA`.

## Details

The trailing dimension of `SrvAge_corr_pars` distinguishes the age AR1
correlation (index 1, used by likelihoods 3 and 4) from the sex
correlation (index 2, used only by likelihood 4). Parameters are
activated conditionally on the fleet's likelihood type and composition
type: aggregated compositions (`type = 0`) activate
`SrvAge_corr_pars_agg` only for likelihood 3; split or joint-by-sex
compositions activate `SrvAge_corr_pars` with index 2 skipped when
`n_sexes = 1`. Fleets using multinomial likelihoods or with no active
age composition data have all parameters mapped to `NA`.

## See also

[`do_SrvAge_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_theta_mapping.md)
for the overdispersion parameter mapping;
[`do_SrvLen_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvLen_corr_pars_mapping.md)
for the analogous length composition correlation mapping.
