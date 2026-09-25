# Set up total and retained fishery selectivity and catchability specifications

Sets the selectivity functional forms, time blocks, continuous time
variation, process error hyperparameters, annual deviations, and the
catchability blocks and estimation structure, for total and retained
selectivity alike. Time variation and blocked selectivity are mutually
exclusive within a fleet. Call after
[`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md).

## Usage

``` r
Setup_Mod_Fishsel_and_Q(
  input_list,
  cont_tv_fish_sel = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ""),
  fish_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ""),
  fish_sel_model,
  Use_fish_q_prior = 0,
  fish_q_prior = NA,
  fish_q_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ""),
  fish_q_type = rep("est", input_list$data$n_fish_fleets),
  fish_q_model = NULL,
  sigma_fish_q_spec = "est_all",
  fish_q_rho_spec = "est_all",
  fish_q_rw_init_sigma = NA,
  fishsel_pe_pars_spec = NULL,
  fish_fixed_sel_pars_spec = NULL,
  fish_q_spec = NULL,
  fish_sel_devs_spec = NULL,
  corr_opt_semipar = NULL,
  Use_fish_selex_prior = 0,
  fish_selex_prior = NULL,
  Use_fish_selex_penalty = 0,
  fish_sel_norm_bins = NULL,
  fish_sel_bin_dev_bins = NULL,
  fishsel_pe_wt = rep(1, input_list$data$n_fish_fleets),
  fishsel_rw_init_sigma = rep(5, input_list$data$n_fish_fleets),
  fishsel_dont_est_dev_first = rep(0, input_list$data$n_fish_fleets),
  cont_tv_fishsel_bin_devs = rep("none", input_list$data$n_fish_fleets),
  fish_selex_penalty = NULL,
  fishsel_devs_shared_bins = NULL,
  fish_selex_type = "age",
  use_fixed_fish_sel = rep(0, input_list$data$n_fish_fleets),
  fish_sel_input = NULL,
  fish_sel_nonpar_est_bins = NULL,
  fish_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
  fish_sel_dbnrml_raw = NULL,
  fish_sel_dbnrml_startbin = NULL,
  cont_tv_ret_sel = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ""),
  ret_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ""),
  ret_sel_model = paste("logist1_Fleet_", 1:input_list$data$n_fish_fleets, sep = ""),
  retsel_pe_pars_spec = NULL,
  ret_fixed_sel_pars_spec = rep("fix_ret_sel_input", input_list$data$n_fish_fleets),
  ret_sel_devs_spec = NULL,
  ret_sel_corr_opt_semipar = NULL,
  Use_ret_selex_prior = 0,
  ret_selex_prior = NULL,
  retsel_devs_shared_bins = NULL,
  retsel_pe_wt = rep(1, input_list$data$n_fish_fleets),
  retsel_rw_init_sigma = rep(5, input_list$data$n_fish_fleets),
  retsel_dont_est_dev_first = rep(0, input_list$data$n_fish_fleets),
  ret_selex_type = "age",
  use_fixed_ret_sel = rep(1, input_list$data$n_fish_fleets),
  ret_sel_input = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages),
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  ret_sel_nonpar_est_bins = NULL,
  ret_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- cont_tv_fish_sel:

  Character vector `[n_fish_fleets]` of continuous time variation per
  fleet, each `"<type>_Fleet_<f>"`: `"none"` (default), `"iid"` or
  `"rw"` on the selectivity parameters, or `"3dmarg"`, `"3dcond"` (3D
  GMRF on the marginal or conditional variance) and `"2dar1"` (separable
  over bin and year). Any fleet other than `"none"` also needs
  `fishsel_pe_pars_spec` and `fish_sel_devs_spec`.

- fish_sel_blocks:

  Character vector of discrete selectivity time blocks per fleet, each
  `"Block_<b>_Year_<s>-<e>_Fleet_<f>"` with `"terminal"` allowed as the
  end year, or `"none_Fleet_<f>"` (default) for one constant block.
  Blocks must not overlap and together must span every model year for
  that fleet. Mutually exclusive with `cont_tv_fish_sel != "none"`.

- fish_sel_model:

  Character vector of the selectivity form per fleet, and optionally per
  block: `"<model>_Fleet_<f>"` or `"<model>_Fleet_<f>_Block_<b>"`. The
  forms are `"logist1"` (\\a\_{50}\\ and slope), `"logist2"`
  (\\a\_{50}\\ and \\a\_{95}\\), `"gamma"` (dome, \\a\_{max}\\ and
  \\\delta\\), `"exponential"` (one power), `"dbnrml"` (double normal,
  six parameters), `"asymplogist1"` and `"asymplogist2"` (the two
  logistics with an asymptote), the three non-parametric forms and
  `"bicubic"`.

  `"nonpar"` is on the logit scale, mean-standardized jointly over years
  and bins so the grand mean of the surface is one. `"nonparlog"` is on
  the log scale, standardized so each year averages to one over
  `*_sel_norm_bins`, leaving only within-year contrasts identified.
  `"nonparfree"` is on the log scale with no standardization,
  \\\exp(\theta)\\, so the values hold the height of the curve as well
  as its shape; this is the form for a data source fit age by age, where
  a free catchability per age and a selectivity estimated at age are one
  quantity written twice, so no catchability is set. Pin one bin, by
  leaving it out of the estimated bins, whenever the mean it multiplies
  is also free.

  `"bicubic"` is a spline over a bin-node by year-node grid, written
  `"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"` with an optional
  `_Block_k`. One form covers a smooth bin by year surface
  (`n_yr_nodes > 1`), a time-invariant bin-only spline
  (`n_yr_nodes == 1`), and a bin-only spline re-fit per block. An
  optional `_SelStyr_<year>` restricts the fit to `SelStyr`:block-end,
  holding earlier years of the block at the `SelStyr` curve, and an
  optional `_NSelBins_<n>` restricts it to the first `n` bins, holding
  the rest at the last fitted bin. See
  [`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md)
  and the model equations vignette.

- Use_fish_q_prior:

  Integer flag, `1` for lognormal priors on catchability. Default `0`.

- fish_q_prior:

  Data frame with columns `region`, `fleet`, `block`, `mu` on the
  natural scale and `sd` on the log scale, one row per
  \\\text{Normal}(\log(\mu), \sigma)\\ prior. Read when
  `Use_fish_q_prior = 1`.

- fish_q_blocks:

  Catchability time blocks per fleet, in the same format as
  `fish_sel_blocks`. Default one constant block.

- fish_q_type:

  Character vector `[n_fish_fleets]` of how catchability is obtained.
  `"est"` (default) estimates `ln_fish_q`, `"arith"` concentrates it out
  as the ratio of mean observed to mean predicted index, and `"geo"`
  does the same on the log scale as `exp(mean(log(obs) - log(pred)))`.
  Both analytic forms use the years with observations only and fix that
  fleet's `ln_fish_q` whatever `fish_q_spec` says. The solve runs within
  each `fish_q_blocks` block, so a blocked catchability gets one solved
  value per block.

- fish_q_model:

  Character vector `[n_fish_fleets]` of the process error on annual
  catchability deviations: `"none"` (default), `"iid"`, `"rw"`, `"ar1"`
  or `"dsem"`, which hands the series to
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md).
  Catchability is then \\\exp(\ln q\_{r,b,f} + \epsilon\_{r,y,f})\\. A
  fleet with deviations cannot also have `fish_q_blocks` or an
  analytically solved `fish_q_type`.

- sigma_fish_q_spec:

  Sharing string for the deviation standard deviation over region and
  fleet: `"est_all"` (default), `"est_shared_r"`, `"est_shared_f"`,
  `"est_shared_r_f"` or `"fix"`.

- fish_q_rho_spec:

  Sharing string for the AR1 correlation, with the same options as
  `sigma_fish_q_spec`. Default `"est_all"`. Only read under
  `fish_q_model = "ar1"`.

- fish_q_rw_init_sigma:

  Standard deviation of the first estimated year of a random walk. `NA`
  (default) starts the walk at zero under its own sigma, which keeps
  `ln_fish_q` as the level of the series.

- fishsel_pe_pars_spec:

  Character vector `[n_fish_fleets]` of the estimation structure for the
  selectivity process error hyperparameters, required when any fleet
  varies continuously. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md).

- fish_fixed_sel_pars_spec:

  Character vector `[n_fish_fleets]` of how the fixed-effect selectivity
  parameters are estimated: `"est_all"`, `"est_shared_r"`,
  `"est_shared_s"`, `"est_shared_r_s"`, `"est_shared_f_x"` or `"fix"`.
  See
  [`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md).

- fish_q_spec:

  Character vector `[n_fish_fleets]` of the catchability estimation
  structure: `"est_all"`, `"est_shared_r"` or `"fix"`. See
  [`do_q_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_q_mapping.md).

- fish_sel_devs_spec:

  Character vector `[n_fish_fleets]` of the estimation structure for the
  annual selectivity deviations, required when any fleet varies
  continuously. See
  [`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md).

- corr_opt_semipar:

  Character vector `[n_fish_fleets]` of which correlation components to
  suppress under 3D GMRF or 2D AR1 time variation. `NA` (default)
  suppresses none, and the cohort options are invalid for `"2dar1"`. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md).

- Use_fish_selex_prior:

  Integer flag, `1` for priors on the selectivity parameters. Default
  `0`.

- fish_selex_prior:

  Data frame with columns `region`, `fleet`, `block`, `sex`, `par`,
  `mu`, `sd` and an optional `type`. `"par"` (the default) is a
  lognormal prior on one fixed selectivity parameter, with `mu` on the
  natural scale and `sd` on the log scale. `"value"` is a normal prior
  on the realized selectivity at one bin, both on the natural scale,
  where `par` names the bin and the value is read at the first model
  year of `block`; that is the ADMB convention of pinning selectivity at
  a reference age near one, which no set of independent parameter priors
  can express. Read when `Use_fish_selex_prior = 1`.

- Use_fish_selex_penalty:

  Integer (0/1). Whether a centering penalty is applied to sets of
  fishery selectivity fixed-effect parameters. Default `0`.

- fish_sel_norm_bins:

  List with one element per fleet naming the bins the mean-one
  standardization averages over, or `NULL` for fleets standardizing over
  every bin. Read under `"nonparlog"` only. A gear whose catchability is
  defined against part of the bin range standardizes over that part, and
  catchability absorbs the difference in scale. Default `NULL`.

- fish_sel_bin_dev_bins:

  List with one element per fleet naming the bins that fleet overrides,
  or `NULL` for none, e.g. `list(1, NULL)`. An overridden bin takes a
  freely estimated annual value \\\exp(\epsilon\_{y,b})\\ in place of
  what the functional form produced, applied after every other
  transformation including standardization, while the rest of the curve
  keeps its parametric shape. Default `NULL`.

- fishsel_pe_wt:

  Numeric vector `[n_fish_fleets]` multiplying the fishery selectivity
  process error likelihood. Default `1`. `0` skips that fleet's process
  error, so the deviations stay estimated but enter the objective only
  through the data and any smoothness or centering penalties. Values
  other than 0 or 1 make an estimated process error sigma
  reinterpretable. Applies to `ln_fishsel_devs` only; the bin-override
  deviations have their own process error.

- fishsel_rw_init_sigma:

  Numeric vector `[n_fish_fleets]` giving the standard deviation of the
  first year of an `"rw"` deviation series. Default `5`, which leaves
  that year effectively free. `NA` instead starts the walk at zero under
  the walk's own estimated sigma, which suits a base curve that already
  describes the first year well.

- fishsel_dont_est_dev_first:

  Integer vector `[n_fish_fleets]` of 0/1, default `0`. Where `1`, that
  fleet's deviations start in year two and the fixed parameters hold
  year one. A non-parametric form has one free base parameter per bin,
  so year one's deviation is that same value written twice with only
  `fishsel_rw_init_sigma` between them, a prior on a level usually meant
  to be free. Refused for the GMRF and 2D AR1 forms, whose deviations
  are a field over years and bins rather than a walk anchored at year
  one.

- cont_tv_fishsel_bin_devs:

  Character vector `[n_fish_fleets]` of the process error on the
  bin-override deviations: `"none"` (default), `"iid"` or `"rw"`. A walk
  has its own estimated sigma per bin.

- fish_selex_penalty:

  Data frame of centering penalties with columns `region`, `fleet`,
  `block`, `sex`, `par` and `wt`, required when
  `Use_fish_selex_penalty = 1`. Each row penalizes
  `wt * (log(mean(exp(pars))))^2` over the parameters named in `par`, a
  single index or a list column of integer vectors. This pins the scalar
  of a non-parametric curve that catchability or fishing mortality would
  otherwise absorb, and is softer than fixing a bin. Meant for parameter
  sets on the log scale. Default `NULL`.

- fishsel_devs_shared_bins:

  List of integer vectors grouping the bins that share one deviation
  series, e.g. `list(1:5, 6:10, 11:30)`. Only read when
  `fish_sel_devs_spec` names an `"est_shared_b"` variant.

- fish_selex_type:

  Character scalar, `"age"` or `"length"`, the bin dim every fishery
  selectivity function is defined over.

- use_fixed_fish_sel:

  Integer vector `[n_fish_fleets]`, `1` to fix fishery selectivity and
  `0` to estimate it.

- fish_sel_input:

  Array of fixed fishery selectivity values
  `[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]`.
  Required when any `use_fixed_fish_sel == 1`.

- fish_sel_nonpar_est_bins:

  Optional bin groupings for non-parametric fishery selectivity,
  structured `[[fleet]][[block]]`, each element a list of integer
  vectors naming the bins that share one estimated parameter. Indices
  are on the bin dim `fish_selex_type` names.

- fish_sel_sex_offset:

  Character vector `[n_fish_fleets]` linking the sexes of a fleet's
  selectivity when `n_sexes > 1`. `"none"` (default) keeps each sex's
  stored parameters its own. `"par"` makes every sex beyond the first
  hold additive offsets on the first sex's stored parameters, so a
  log-scale parameter's natural value is the first sex's times
  \\e^{\delta}\\; offsets fixed at zero reproduce sex-shared parameters.
  `"scale"` keeps each sex's own parameters and adds a constant
  log-scale offset on the whole realized curve,
  `exp(ln_fishsel_sex_scale)`, per region, block and sex, which may
  exceed one and is refused for the non-parametric forms and
  semi-parametric time variation, whose standardization would cancel it.
  `"apical"` has the double normal build its limbs up to
  `exp(ln_*sel_sex_scale)` rather than one, so the offset moves the
  middle of the curve and leaves its ends where that sex's own
  parameters put them. `"par_apical"` and `"par_scale"` combine a par
  offset with each.

- fish_sel_dbnrml_raw:

  `NULL` (default) or a 0/1 matrix `[n_fish_fleets x 2]` for fleets on
  the double normal: column one leaves the ascending limb a raw Gaussian
  instead of anchoring it to `p5` at the first bin, column two does the
  same for the descending limb and `p6`.

- fish_sel_dbnrml_startbin:

  `NULL` (default) or an integer vector `[n_fish_fleets]`, the bin each
  fleet's double normal anchors its ascending limb at. Bins below it
  take the squared ratio of their bin to it times the selectivity there,
  which is Stock Synthesis's convention when the compositions start
  above the population's first length bin.

- cont_tv_ret_sel:

  Continuous time variation on retention, with the options and
  requirements of `cont_tv_fish_sel`. Any fleet other than `"none"` also
  needs `retsel_pe_pars_spec` and `ret_sel_devs_spec`.

- ret_sel_blocks:

  Discrete retention time blocks per fleet, in the format of
  `fish_sel_blocks` and mutually exclusive with
  `cont_tv_ret_sel != "none"`.

- ret_sel_model:

  Retention selectivity form per fleet and block, with the same syntax
  and the same set of forms as `fish_sel_model`.

- retsel_pe_pars_spec:

  Estimation structure for the retention process error hyperparameters,
  as `fishsel_pe_pars_spec`.

- ret_fixed_sel_pars_spec:

  How the retention fixed-effect parameters are estimated, with the
  options of `fish_fixed_sel_pars_spec`.

- ret_sel_devs_spec:

  Estimation structure for the annual retention deviations, as
  `fish_sel_devs_spec`.

- ret_sel_corr_opt_semipar:

  Which correlation components to suppress under semi-parametric
  retention time variation, as `corr_opt_semipar`.

- Use_ret_selex_prior:

  Integer flag, `1` for priors on the retention selectivity parameters.
  Default `0`.

- ret_selex_prior:

  Data frame with the columns and the optional `type` of
  `fish_selex_prior`.

- retsel_devs_shared_bins:

  Bins sharing one retention deviation series, as
  `fishsel_devs_shared_bins`.

- retsel_pe_wt:

  Per-fleet multiplier on the retention process error likelihood, as
  `fishsel_pe_wt`. Default `1`.

- retsel_rw_init_sigma:

  Standard deviation of the first year of an `"rw"` retention deviation
  series, as `fishsel_rw_init_sigma`. Default `5`.

- retsel_dont_est_dev_first:

  Whether each fleet's retention deviations start in year two, as
  `fishsel_dont_est_dev_first`. Default `0`.

- ret_selex_type:

  Character scalar, `"age"` or `"length"`, the bin dim every retention
  selectivity function is defined over.

- use_fixed_ret_sel:

  Integer vector `[n_fish_fleets]`, `1` to fix retention selectivity and
  `0` to estimate it.

- ret_sel_input:

  Array of fixed retention values
  `[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]`.

- ret_sel_nonpar_est_bins:

  Optional bin groupings for non-parametric retention, structured
  `[[fleet]][[block]]`, each element a list of bin index vectors
  defining grouped parameters.

- ret_sel_sex_offset:

  Character vector `[n_fish_fleets]` linking the sexes of a fleet's
  retention curve, with the options of `fish_sel_sex_offset`. Default
  `"none"`. Retention is a fraction, so a scale offset only makes sense
  where the scaled curve stays at or below one.

- ...:

  Optional starting values for the selectivity parameters.

## Value

`input_list` with `$data`, `$par` and `$map` updated: the parsed integer
arrays for `cont_tv_fish_sel`, `fish_sel_blocks`, `fish_sel_model` and
`fish_q_blocks`, the starting values for all four parameter groups, and
their factor maps.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
