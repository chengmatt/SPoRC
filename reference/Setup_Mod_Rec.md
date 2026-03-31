# Setup model objects for specifying recruitment module and associated processes

Setup model objects for specifying recruitment module and associated
processes

## Usage

``` r
Setup_Mod_Rec(
  input_list,
  rec_model,
  rec_dd = NULL,
  rec_lag = 1,
  Use_h_prior = 0,
  h_prior = NULL,
  Use_Rec_prop_Prior = 0,
  Rec_prop_prior = NULL,
  do_rec_bias_ramp = 0,
  bias_year = NA,
  max_bias_ramp_fct = 1,
  sigmaR_switch = 1,
  dont_est_recdev_last = 0,
  init_age_strc = 2,
  equil_init_age_strc = 1,
  init_F_prop = 0,
  sigmaR_spec = NULL,
  InitDevs_spec = NULL,
  RecDevs_spec = NULL,
  h_spec = NULL,
  t_spawn = 0,
  sexratio_spec = "fix",
  sexratio_blocks = c(paste("none_Region_", c(1:input_list$data$n_regions), sep = "")),
  ...
)
```

## Arguments

- input_list:

  List containing data, parameters, and map lists used by the model.

- rec_model:

  Character string specifying the recruitment model. Options are:

  - `"mean_rec"`: Recruitment is a fixed mean value.

  - `"bh_rec"`: Beverton-Holt recruitment with steepness parameter.

- rec_dd:

  Character string specifying recruitment density dependence, options:
  `"local"`, `"global"`, or `NULL`.

- rec_lag:

  Integer specifying the recruitment lag duration relative to spawning
  stock biomass (SSB).

- Use_h_prior:

  Integer flag (0 or 1) indicating whether to apply a prior on steepness
  `h`.

- h_prior:

  Data frame specifying beta prior distributions for the \`h_trans\`
  parameters. Must include the following columns: - \`region\`: Integer
  region index corresponding to the element in \`h_trans\` being
  penalized. - \`mu\`: Mean of the prior in normal space (used to
  calculate the corresponding beta distribution). - \`sd\`: Standard
  deviation of the prior in normal space. For each row, a beta
  distribution is scaled to the interval \[0.2, 1\], and the
  corresponding element of \`h_trans\` is transformed to that scale and
  penalized using the log-density from the beta distribution.

- Use_Rec_prop_Prior:

  Integer flag (0 or 1) indicating whether to apply a prior on
  recruitment proportions.

- Rec_prop_prior:

  Scalar or array specifying prior values for recruitment proportion
  parameters. If scalar, a constant uniform prior is applied across all
  dimensions.

- do_rec_bias_ramp:

  Integer flag (0 or 1) indicating whether to apply a recruitment bias
  correction ramp.

- bias_year:

  Numeric vector of length 4 defining the recruitment bias ramp periods:

  - Element 1: End year of no bias correction period.

  - Element 2: End year of ascending bias ramp period.

  - Element 3: End year of full bias correction period.

  - Element 4: Start year of final no bias correction period.

  For example, with 65 years total, `c(21, 31, 60, 64)` means:

  - Years 1–21: No bias correction.

  - Years 22–31: Ascending bias correction.

  - Years 32–60: Full bias correction.

  - Years 61–63: Descending bias ramp.

  - Years 64–65: No bias correction.

- max_bias_ramp_fct:

  Numeric specifying the maximum bias correction to apply to the
  recruitment bias ramp (should be between 0 and 1)

- sigmaR_switch:

  Integer year indicating when `sigmaR` switches from early to late
  values (0 disables switching).

- dont_est_recdev_last:

  Integer specifying how many of the most recent recruitment deviations
  to not estimate. Default is 0.

- init_age_strc:

  Integer flag specifying initialization of initial age structure:

  - `0`: Initialize by iteration.

  - `1`: Initialize using a scalar geometric series w/o any movement in
    all groups (does not account for movement in all group).

  - `2`: Initialize using a matrix geometric series (accounts for
    movement; default).

  - `3`: Initialize using a scalar geometric series w/o any movement in
    only plus groups (does not account for movement in plus group).

- equil_init_age_strc:

  Integer flag specifying how initial age structure deviations should be
  initialized. Default is stochastic for all ages except the recruitment
  age and the plus group.

  - `0`: Equilibrium initial age structure.

  - `1`: Stochastic initial age structure for all ages, except for the
    plus group, which follows equilibrium calculations (geometric
    series)

  - `2`: Stochastic initial age structure for all ages

- init_F_prop:

  Numeric value specifying the initial fishing mortality proportion
  relative to mean fishing mortality for initializing age structure.

- sigmaR_spec:

  Character string specifying estimation of recruitment variability
  (`sigmaR`):

  - `NULL or "est_all"`: Estimate separate `sigmaR` for early and late
    periods.

  - `"est_shared"`: Estimate one `sigmaR` shared across periods.

  - `"fix"`: Fix both `sigmaR` values.

  - `"fix_early_est_late"`: Fix early `sigmaR`, estimate late `sigmaR`.

- InitDevs_spec:

  Character string specifying estimation of initial age deviations:

  - `NULL`: Estimate deviations for all ages and regions.

  - `"est_shared_r"`: Estimate deviations shared across regions.

  - `"fix"`: Fix all deviations.

- RecDevs_spec:

  Character string specifying recruitment deviation estimation:

  - `NULL`: Estimate deviations for all regions and years.

  - `"est_shared_r"`: Estimate deviations shared across regions (global
    recruitment deviations).

  - `"fix"`: Fix all recruitment deviations.

- h_spec:

  Character string specifying steepness estimation:

  - `NULL`: Estimate steepness for all regions if
    `rec_model == "bh_rec"`.

  - `"est_shared_r"`: Estimate steepness shared across regions.

  - `"fix"`: Fix steepness values.

  If `rec_model == "mean_rec"`, steepness is fixed.

- t_spawn:

  Numeric fraction specifying spawning timing within the year.

- sexratio_spec:

  Character string specifying sex ratio estimation scheme:

  - `"est_all"` estimates sex ratio for all blocks and regions
    independently

  - `"est_shared_r"` estimates sex ratio shared across regions but
    varying by block

  - `"fix"` fixes all sex ratio (no estimation)

- sexratio_blocks:

  Character vector specifying blocks of years and regions for sex ratio.
  Format examples:

  - `"Block_1_Year_1-15_Region_1"`

  - `"Block_2_Year_16-terminal_Region_2"`

  - `"none_Region_3"` (means no block, constant for that region; this is
    the default option)

- ...:

  Additional arguments specifying starting values for recruitment
  parameters such as `ln_global_R0`, `Rec_prop`, `h`, `ln_InitDevs`,
  `ln_RecDevs`, and `ln_sigmaR`.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
