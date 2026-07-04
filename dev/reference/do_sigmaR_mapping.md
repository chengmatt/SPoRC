# Map recruitment variability (sigma_R) parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `ln_sigmaR`, the log-scale
standard deviation of recruitment deviations. The `ln_sigmaR` array has
dimensions `[2 x n_pop x n_regions]`, where the first index
distinguishes initial age-structure deviations (`i = 1`) from annual
recruitment deviations (`i = 2`).

## Usage

``` r
do_sigmaR_mapping(input_list, sigmaR_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions.

- sigmaR_spec:

  Character string specifying the estimation structure for `ln_sigmaR`.
  One of:

  `"est_all"`

  :   Estimate `ln_sigmaR` separately for both the initial deviation
      period (`i = 1`) and the annual deviation period (`i = 2`),
      respecting the density-dependence sharing rules described above.

  `"est_shared_all"`

  :   Single `ln_sigmaR` shared across both deviation periods, all
      populations, and all regions.

  `"fix_early_est_late"`

  :   Fix `ln_sigmaR` for the initial deviation period (`i = 1`) at its
      starting value; estimate `ln_sigmaR` for the annual deviation
      period (`i = 2`). Useful when initial age-structure uncertainty is
      assumed known.

  `"fix"`

  :   Fix all `ln_sigmaR` parameters at their starting values (all
      mapped to `NA`).

## Value

The input `input_list` with `$map$ln_sigmaR` set to a factor vector of
length `prod(dim(par$ln_sigmaR))`. Active parameters receive sequential
integer indices; fixed parameters are `NA`.

## Details

Sharing behaviour adapts to the density-dependence structure: under
local density dependence (`rec_dd = 0`) with a single population, each
region receives its own `ln_sigmaR` estimate; under global density
dependence or with multiple populations, a single `ln_sigmaR` is shared
across regions within each population.
