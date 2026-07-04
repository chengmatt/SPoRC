# Map fishery length composition overdispersion parameters

Analogous to
[`do_FishAge_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_FishAge_theta_mapping.md)
but for length compositions. Constructs factor maps for
`ln_FishLen_theta` and `ln_FishLen_theta_agg` based on
`FishLenComps_Type`, `FishLenComps_LikeType`, and `UseFishLenComps`.
Parameters are mapped to `NA` for fleets using multinomial likelihoods
or with no observed length compositions.

## Usage

``` r
do_FishLen_theta_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `FishLenComps_Type`, `FishLenComps_LikeType`, and `UseFishLenComps` to
  be set by
  [`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md).

## Value

The input `input_list` with `$map$ln_FishLen_theta` and
`$map$ln_FishLen_theta_agg` set to factor vectors.
