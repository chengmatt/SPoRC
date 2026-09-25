# Map recruitment variability (sigma_R) parameters

Builds the factor map for `ln_sigmaR` `[2 x n_pop x n_regions]`, where
the first index is the initial deviation period and the second the
annual one. The recruitment penalty reads each region's own slot, and a
population's non-natal region slots are mapped off under
`rec_region_prop_spec = 1` with several populations. Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

## Usage

``` r
do_sigmaR_mapping(input_list, sigmaR_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- sigmaR_spec:

  Estimation structure for `ln_sigmaR`. `"est_all"` estimates each
  period, population and region separately. `"est_shared_r"` estimates
  per period and population, shared across regions. `"est_shared_all"`
  gives one value for everything. `"fix_early_est_late"` holds the
  initial period at its starting value and estimates the annual one per
  population and region. `"fix"` holds every value.

## Value

`input_list` with `$map$ln_sigmaR` set to a factor vector of length
`prod(dim(par$ln_sigmaR))`. Active parameters take sequential integers,
fixed ones are `NA`.
