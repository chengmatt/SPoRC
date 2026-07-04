# Map population survey length composition correlation parameters

Constructs factor maps for `SrvLen_pop_corr_pars` (region- and
sex-specific AR1 and sex correlation parameters) and
`SrvLen_pop_corr_pars_agg` (aggregated correlation parameters) for 1D
and 2D logistic-normal len composition likelihoods.

## Usage

``` r
do_SrvLen_pop_corr_pars_mapping(input_list)
```

## Arguments

- input_list:

  Named list containing `data`, `par`, and `map` components.

## Value

The input `input_list` with elements `map\$SrvLen_pop_corr_pars` and
`map\$SrvLen_pop_corr_pars_agg` set to factor vectors.

## Details

Parameters are activated only when `SrvLenComps_pop_LikeType` is in
`c(3, 4)`. These correspond to the 1D and 2D logistic-normal
likelihoods. All other likelihoods map correlation parameters to `NA`.

For the 2D logistic-normal (`LikeType == 4`), both trailing elements of
the `[,,,,2]` slice are activated: element 1 for the len AR1 coefficient
and element 2 for the sex correlation (skipped when `n_sexes == 1`).
