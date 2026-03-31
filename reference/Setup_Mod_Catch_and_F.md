# Setup fishing mortality and catch observations

Setup fishing mortality and catch observations

## Usage

``` r
Setup_Mod_Catch_and_F(
  input_list,
  ObsCatch,
  catch_units = array("biom", dim = c(input_list$data$n_regions,
    input_list$data$n_fish_fleets)),
  Catch_Type = array(1, dim = c(length(input_list$data$years),
    input_list$data$n_fish_fleets)),
  UseCatch,
  Use_F_pen = 1,
  est_all_regional_F = 1,
  sigmaC_spec = "fix",
  sigmaC_agg_spec = "fix",
  sigmaF_spec = "fix",
  sigmaF_agg_spec = "fix",
  ...
)
```

## Arguments

- input_list:

  A list containing data, parameters, and map lists used by the model.

- ObsCatch:

  Numeric array of observed catches, dimensioned
  `[n_regions, n_years, n_fish_fleets]`.

- catch_units:

  Catch units - Array dimensioned by n_regions x n_fish_fleets

  - `"abd"`: Catch units in abundance

  - `"biom"`: Catch units in biomass (default)

- Catch_Type:

  Integer matrix with dimensions `[n_years, n_fish_fleets]`, specifying
  catch data types:

  - `0`: Use aggregated catch data for the year.

  - `1`: Use region-specific catch data for the year.

- UseCatch:

  Indicator array `[n_regions, n_years, n_fish_fleets]` specifying
  whether to include catch data in the fit:

  - `0`: Do not use catch data.

  - `1`: Use catch data and fit.

- Use_F_pen:

  Integer flag indicating whether to apply a fishing mortality penalty:

  - `0`: Do not apply penalty.

  - `1`: Apply penalty.

- est_all_regional_F:

  Integer flag indicating whether all regional fishing mortality
  deviations are estimated:

  - `0`: Some fishing mortality deviations are aggregated across
    regions.

  - `1`: All fishing mortality deviations are regional.

- sigmaC_spec:

  Character string specifying observation error structure for catch
  data. Default behavior fixes `sigmaC` at a starting value of `1e-3`
  (log-scale `ln_sigmaC = log(1e-3)`) for all regions, years, and
  fleets. Other options include:

  - `"est_shared_f"`: Estimate `sigmaC` shared across fishery fleets,
    unique by region and year.

  - `"est_shared_r"`: Estimate `sigmaC` shared across regions, unique by
    fleet and year.

  - `"est_shared_y"`: Estimate `sigmaC` shared across years, unique by
    region and fleet.

  - `"est_shared_r_f"`: Estimate `sigmaC` shared across regions and
    fleets, unique by year.

  - `"est_shared_f_y"`: Estimate `sigmaC` shared across fleets and
    years, unique by region.

  - `"est_shared_r_y"`: Estimate `sigmaC` shared across regions and
    years, unique by fleet.

  - `"est_shared_r_y_f"`: Estimate single `sigmaC` shared across
    regions, years, and fleets.

  - `"fix"`: Fix `sigmaC` at the starting value.

  - `"est_all"`: Estimate separate `sigmaC` for each region, year, and
    fleet combination.

- sigmaC_agg_spec:

  Character string specifying process error structure for aggregated
  catch observation error. Default fixes `sigmaC_agg` at the starting
  value (log-scale `sigmaC_agg`). Other options include:

  - `"est_shared_f"`: Estimate `sigmaC_agg` shared across fishery
    fleets, unique by year.

  - `"est_shared_y"`: Estimate `sigmaC_agg` shared across years, unique
    by fleet.

  - `"est_shared_y_f"`: Estimate single `sigmaC_agg` shared across years
    and fleets.

  - `"fix"`: Fix at the starting value.

  - `"est_all"`: Estimate separate parameters for each year and fleet
    combination.

- sigmaF_spec:

  Character string specifying process error structure for fishing
  mortality. Default fixes `sigmaF` at `1` on the log scale (i.e.,
  `ln_sigmaF = 0`). Other options include:

  - `"est_shared_f"`: Estimate `sigmaF` shared across fishery fleets.

  - `"est_shared_r"`: Estimate `sigmaF` shared across regions but unique
    by fleet.

  - `"est_shared_r_f"`: Estimate `sigmaF` shared across regions and
    fleets.

  - `"fix"`: Fix `sigmaF` at the starting value.

  - `"est_all"`: Estimate separate `sigmaF` for each region and fleet.

- sigmaF_agg_spec:

  Character string specifying process error structure for aggregated
  fishing mortality. Default fixes `sigmaF_agg` at the starting value
  (log-scale `ln_sigmaF_agg`). Other options include:

  - `"est_shared_f"`: Estimate `sigmaF_agg` shared across fishery
    fleets.

  - `"fix"`: Fix at the starting value.

  - `"est_all"`: Estimate separate parameters for each fishery fleet.

- ...:

  Additional arguments specifying starting values for `ln_sigmaC`,
  `ln_sigmaF`, and `ln_sigmaF_agg`.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
