# Map fishery length composition correlation parameters

Analogous to
[`do_FishAge_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_FishAge_corr_pars_mapping.md)
but for length compositions. Constructs factor maps for
`FishLen_corr_pars` and `FishLen_corr_pars_agg` for 1D and 2D
logistic-normal length composition likelihoods (`FishLenComps_LikeType`
is in `c(3, 4)`).

## Usage

``` r
do_FishLen_corr_pars_mapping(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

## Value

The input `input_list` with `$map$FishLen_corr_pars` and
`$map$FishLen_corr_pars_agg` set to factor vectors.
