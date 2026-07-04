# Map discarded fishery length composition overdispersion parameters

Constructs factor maps for `ln_FishLen_discard_theta` (fleet- region-
and sex-specific overdispersion) and `ln_FishLen_discard_theta_agg`
(aggregated overdispersion) based on the composition type and likelihood
specified in `$data$FishLenComps_discard_Type` and
`$data$FishLenComps_discard_LikeType`. Parameters are mapped to `NA` for
fleets using multinomial likelihoods (`LikeType == 0`) or with no
observed len compositions.

## Usage

``` r
do_FishLen_discard_theta_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `FishLenComps_discard_Type`, `FishLenComps_discard_LikeType`, and
  `UseFishLenComps_discard` to be set by
  [`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md).

## Value

The input `input_list` with `$map$ln_FishLen_discard_theta` and
`$map$ln_FishLen_discard_theta_agg` set to factor vectors. Active
parameters receive sequential integer indices; inactive parameters are
`NA`.
