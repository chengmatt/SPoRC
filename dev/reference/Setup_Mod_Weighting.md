# Set likelihood and penalty weights for the estimation model

Assigns the \\\lambda\\ multiplier on each likelihood component and
penalty term, which is how a noisy data source is down-weighted, how an
iterative reweighting such as Francis is applied, and how a component is
switched off entirely, with a weight of `0`. Call after every data setup
function.

## Usage

``` r
Setup_Mod_Weighting(
  input_list,
  addtocomp = 0.001,
  comp_const_obs = 1,
  addtofishidx = 1e-04,
  addtosrvidx = 1e-04,
  addtotag = 1e-10,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1,
  Wt_Catch_pop = 1,
  Wt_FishIdx_pop = 1,
  Wt_SrvIdx_pop = 1,
  Wt_Rec = 1,
  Wt_Init_Rec = NULL,
  Wt_F = 1,
  Wt_Tagging = 1,
  Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_FishAgeComps_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_FishLenComps_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_SrvLenComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_Discard = 1,
  Wt_Discard_pop = 1,
  Wt_D = 1,
  Wt_FishAgeComps_discard = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_FishLenComps_discard = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_FishAgeComps_discard_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_FishLenComps_discard_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_Fish_caal = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens),
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_Srv_caal = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens),
    input_list$data$n_sexes, input_list$data$n_srv_fleets)),
  fish_sel_pen_wts = NULL,
  ret_sel_pen_wts = NULL,
  srv_sel_pen_wts = NULL
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- addtocomp:

  Small constant added to the composition proportions to avoid `log(0)`.
  Default `1e-3`. Ignored by the logistic normal, which handles zeros
  itself.

- comp_const_obs:

  Integer switch for where `addtocomp` enters the multinomial, not a
  constant to tune. `1` (default) adds it to the observed proportions
  that weight the likelihood as well as inside the logarithms, so the
  likelihood is stationary at `pred = obs`; `0` weights by the raw
  observed proportions. With a Dirichlet-multinomial conditional
  age-at-length fleet, `1` warns, since the constant biases theta upward
  when most age bins in a length bin are structurally empty.

- addtofishidx, addtosrvidx:

  Small constants added to the fishery and survey indices. Default
  `1e-4`.

- addtotag:

  Small constant added to the tag recovery observations. Default
  `1e-10`.

- Wt_Catch, Wt_FishIdx:

  Weights on the catch and fishery index likelihoods, a scalar or an
  array `[n_regions × n_years × n_seas × n_fish_fleets]`. Default `1`.

- Wt_SrvIdx:

  Weight on the survey index likelihood, a scalar or an array
  `[n_regions × n_years × n_seas × n_srv_fleets]`. Default `1`.

- Wt_Catch_pop, Wt_FishIdx_pop:

  The population-specific catch and fishery index weights, a scalar or
  an array `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.
  Default `1`.

- Wt_SrvIdx_pop:

  The population-specific survey index weight, a scalar or an array
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`. Default `1`.

- Wt_Rec:

  Weight on the recruitment deviation penalty, a scalar or an array
  `[n_pop × n_regions × n_est_rec_devs]`, where the third dim is
  `ln_RecDevs`'s own rather than the number of years, since
  `dont_est_recdev_last` and `n_proj_yrs_devs` both move it. Default
  `1`. A zero leaves a deviation estimated but takes it out of the
  penalty, which is how a stock-recruit relationship is fit over a
  window of years while recruitment stays free in every year;
  `dont_est_recdev_last` instead removes the deviations, so recruitment
  reverts to the deterministic prediction.

- Wt_Init_Rec:

  Weight on the initial age deviation penalty, a scalar or an array
  `[n_pop × n_regions × (n_ages - 1) × n_sexes]`. `NULL` (default) takes
  `Wt_Rec` when that is a scalar; supply it explicitly when `Wt_Rec` is
  an array, since the two penalties are dimensioned differently.

- Wt_F:

  Scalar weight on the fishing mortality deviation penalty. Default `1`.

- Wt_Tagging:

  Scalar weight on the tag recovery likelihood. Default `1`.

- Wt_FishAgeComps, Wt_FishLenComps, Wt_FishAgeComps_discard,
  Wt_FishLenComps_discard:

  Weights on the fishery and discard composition likelihoods, a scalar
  or an array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. Default
  one everywhere.

- Wt_SrvAgeComps, Wt_SrvLenComps:

  Weights on the survey composition likelihoods, a scalar or an array
  `[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`. Default one
  everywhere.

- Wt_FishAgeComps_pop, Wt_FishLenComps_pop, Wt_FishAgeComps_discard_pop,
  Wt_FishLenComps_discard_pop:

  The population-specific fishery and discard composition weights, a
  scalar or an array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  Default one everywhere.

- Wt_SrvAgeComps_pop, Wt_SrvLenComps_pop:

  The population-specific survey composition weights, a scalar or an
  array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`.
  Default one everywhere.

- Wt_Discard:

  Weight on the aggregated discard amount or fraction likelihood, a
  scalar or an array `[n_regions × n_years × n_seas × n_fish_fleets]`.
  Default `1`.

- Wt_Discard_pop:

  The population-specific discard weight, a scalar or an array with a
  leading `n_pop` dim. Default `1`.

- Wt_D:

  Scalar weight on the discard mortality rate deviation penalty. Default
  `1`.

- Wt_Fish_caal:

  Weight on the fishery conditional age-at-length likelihood,
  multiplying each length bin's input sample size. Array
  `[n_regions x n_years x n_seas x n_lens x n_sexes x n_fish_fleets]`,
  the shape of `ISS_Fish_caal`. Default one everywhere.

- Wt_Srv_caal:

  The survey counterpart, with `n_srv_fleets` last. Default one
  everywhere.

- fish_sel_pen_wts:

  `NULL` (default), or a named numeric vector or list weighting any
  subset of six selectivity smoothness penalties, which are evaluated on
  the fleet's realized selectivity by bin and year surface and so apply
  to any functional form: `"smooth_bin_curve"` and `"smooth_bin_diff"`
  are the second and first difference across bins, `"smooth_yr_diff"`
  and `"smooth_yr_curve"` the same across years, `"smooth_dome"`
  penalizes non-monotonicity across bins, and `"smooth_mean_center"`
  regularizes each year's mean. See
  [`resolve_sel_pen_wts`](https://chengmatt.github.io/SPoRC/dev/reference/resolve_sel_pen_wts.md)
  and
  [`Get_Selex_Smoothness_Penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex_Smoothness_Penalty.md).
  A name left out is `0`. Each weight may instead be a vector with one
  value per model year, so a penalty can act in some years only or at a
  different strength in each. The specification may also hold
  `"bin_range"`, the first and last bin the penalties act over. Pass an
  unnamed list of per-fleet specifications to give each fleet its own.
  Call after `Setup_Mod_Fishsel_and_Q`.

- ret_sel_pen_wts:

  As `fish_sel_pen_wts`, for retained fishery selectivity.

- srv_sel_pen_wts:

  As `fish_sel_pen_wts`, for survey selectivity. Call after
  `Setup_Mod_Srvsel_and_Q`.

## Value

`input_list` with every weight stored in `$data` under its own name.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
