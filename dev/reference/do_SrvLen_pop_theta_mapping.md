# Map population survey length composition overdispersion parameters

Constructs factor maps for `ln_SrvLen_pop_theta` (fleet- region- and
sex-specific overdispersion) and `ln_SrvLen_pop_theta_agg` (aggregated
overdispersion) based on the composition type and likelihood specified
in `$data$SrvLenComps_pop_Type` and `$data$SrvLenComps_pop_LikeType`.
Parameters are mapped to `NA` for fleets using multinomial likelihoods
(`LikeType == 0`) or with no observed len compositions.

## Usage

``` r
do_SrvLen_pop_theta_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `SrvLenComps_pop_Type`, `SrvLenComps_pop_LikeType`, and
  `UseSrvLenComps_pop` to be set by
  [`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md).

## Value

The input `input_list` with `$map$ln_SrvLen_pop_theta` and
`$map$ln_SrvLen_pop_theta_agg` set to factor vectors. Active parameters
receive sequential integer indices; inactive parameters are `NA`.
