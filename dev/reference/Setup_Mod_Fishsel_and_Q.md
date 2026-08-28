# Set up total and retained fishery selectivity and catchability specifications

Configures all aspects of fishery selectivity and catchability for the
estimation model: functional forms, time blocks, continuous time-varying
structures, process error hyperparameters, annual deviations, and
catchability blocks and estimation structure. Must be called after
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
  fish_q_cov_dat = NULL,
  fish_q_formula = NULL,
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

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists.

- cont_tv_fish_sel:

  Character vector of length `n_fish_fleets` specifying continuous
  time-varying selectivity per fleet. Each element must be
  `"<type>_Fleet_<f>"`. Valid types:

  `"none"`

  :   No continuous time-variation (default).

  `"iid"`

  :   IID annual deviations on selectivity parameters.

  `"rw"`

  :   Random walk in selectivity parameters over time.

  `"3dmarg"`

  :   3D GMRF with marginal variance parameterization.

  `"3dcond"`

  :   3D GMRF with conditional variance parameterization.

  `"2dar1"`

  :   2D separable AR1 in bin and year dimensions.

  If any fleet has `cont_tv_fish_sel != "none"`, both
  `fishsel_pe_pars_spec` and `fish_sel_devs_spec` must also be provided.

- fish_sel_blocks:

  Character vector defining discrete selectivity time blocks per fleet.
  Each element follows `"Block_<b>_Year_<s>-<e>_Fleet_<f>"` or
  `"Block_<b>_Year_<s>-terminal_Fleet_<f>"`. Use `"none_Fleet_<f>"`
  (default) for a single constant block. Blocks must be non-overlapping
  and together span all model years for the specified fleet. Mutually
  exclusive with `cont_tv_fish_sel != "none"` for the same fleet.

- fish_sel_model:

  Character vector specifying the selectivity functional form for each
  fleet (and optionally each time block). Each element must follow one
  of:

  - `"<model>_Fleet_<f>"`: single form for all years of fleet `f`.

  - `"<model>_Fleet_<f>_Block_<b>"`: form specific to block `b` of fleet
    `f`, as defined in `fish_sel_blocks`.

  Available models:

  `"logist1"`

  :   Logistic with \\a\_{50}\\ and slope \\k\\ (2 parameters).

  `"logist2"`

  :   Logistic with \\a\_{50}\\ and \\a\_{95}\\ (2 parameters).

  `"gamma"`

  :   Dome-shaped gamma with \\a\_{max}\\ and \\\delta\\ (2 parameters).

  `"exponential"`

  :   Exponential with a single power parameter (1 parameter).

  `"dbnrml"`

  :   Double-normal with 6 parameters.

  `"nonpar"`

  :   Non-parametric over discrete age or length bins, on the logit
      scale, then mean-standardized jointly over years and bins so the
      grand mean of the surface is one. Bins may be grouped through the
      non-parametric bin mapping. No fixed functional form is imposed.

  `"nonparlog"`

  :   Non-parametric on the log scale, standardized so each year's
      selectivity averages to one over `*_sel_norm_bins`. Only
      within-year contrasts are identified; the level is absorbed by
      catchability or fishing mortality.

  `"nonparfree"`

  :   Non-parametric on the log scale with no standardization,
      \\\exp(\theta)\\, so the values carry the height of the curve as
      well as its shape. This is the form for a stream fit age by age: a
      free catchability per age and a selectivity estimated at age are
      one quantity written two ways, so the whole age multiplier lives
      here and no catchability is set. Pin one bin, by leaving it out of
      the estimated bins, whenever the mean it multiplies is also free.

  `"asymplogist1"`

  :   Logistic selectivity with \\a\_{50}\\ and slope \\k\\ and
      asymptotic control (3 parameters).

  `"asymplogist2"`

  :   Logistic selectivity with with \\a\_{50}\\ and \\a\_{95}\\ and
      asymptotic control (3 parameters).

  `"bicubic"`

  :   Bicubic spline over a bin-node x year-node grid (see
      [`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md),
      `Selex_Model == 8`). Specified as
      `"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"` (optionally
      with `_Block_k`). One generalized form covers a smooth bin x year
      surface (`n_yr_nodes > 1`), a time-invariant bin-only spline
      (`n_yr_nodes == 1`), or a bin-only spline re-fit independently per
      year-block (`n_yr_nodes == 1` within each of several blocks
      defined via `fish_sel_blocks`). An optional `_SelStyr_<year>`
      suffix (a calendar year within the block) restricts the actual
      spline fit to `SelStyr`:block-end; years within the block before
      `SelStyr` are held constant at the `SelStyr` year's fitted curve,
      rather than fitting the surface over the whole block. An optional
      `_NSelBins_<n>` suffix restricts the actual spline fit to the
      first `n` bins (ages or lengths, per `fish_selex_type`); bins
      beyond `n` are held constant at the last fitted bin's curve.

  See the model equations vignette for mathematical definitions.

- Use_fish_q_prior:

  Integer flag. `1` = apply lognormal priors to catchability; `0` = no
  priors (default). Requires `fish_q_prior`.

- fish_q_prior:

  Data frame of catchability prior hyperparameters. Required columns:
  `region`, `fleet`, `block` (block index), `mu` (prior mean on natural
  scale), `sd` (prior SD on log scale). Each row specifies a
  \\\text{Normal}(\log(\mu), \sigma)\\ prior for one catchability
  parameter. Only used when `Use_fish_q_prior = 1`.

- fish_q_blocks:

  Character vector defining catchability time blocks per fleet, using
  the same format as `fish_sel_blocks`. Default `"none_Fleet_<f>"` gives
  a single constant block.

- fish_q_type:

  Character vector of length `n_fish_fleets` controlling how
  catchability is obtained. `"est"` (default) estimates `ln_fish_q`.
  `"arith"` concentrates it out of the likelihood as the ratio of mean
  observed to mean predicted index, and `"geo"` does the same on the log
  scale as `exp(mean(log(obs) - log(pred)))`. Both analytic forms solve
  one catchability per region and fleet using only the years with
  observations, ignore any block structure, and fix that fleet's
  `ln_fish_q` regardless of `fish_q_spec`.

- fish_q_cov_dat:

  Named list of numeric vectors (length = `n_years`) containing the
  covariate time series referenced in `fish_q_formula`. All vectors must
  be the same length and contain no missing values; set values to `0`
  for years when the fishery index is not active. Default `NULL`.

- fish_q_formula:

  Named list of one-sided formulas, one element per fleet requiring
  catchability covariates, referencing series in `fish_q_cov_dat`.
  `NULL` (default) excludes covariate effects.

- fishsel_pe_pars_spec:

  Character vector of length `n_fish_fleets` specifying the estimation
  structure for selectivity process error hyperparameters. Required when
  any fleet has continuous time-variation. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md)
  for all options.

- fish_fixed_sel_pars_spec:

  Character vector of length `n_fish_fleets` specifying how fixed-effect
  selectivity parameters are estimated. See
  [`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md)
  for all options (`"est_all"`, `"est_shared_r"`, `"est_shared_s"`,
  `"est_shared_r_s"`, `"est_shared_f_x"`, `"fix"`).

- fish_q_spec:

  Character vector of length `n_fish_fleets` specifying catchability
  estimation structure. See
  [`do_q_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_q_mapping.md)
  for options (`"est_all"`, `"est_shared_r"`, `"fix"`).

- fish_sel_devs_spec:

  Character vector of length `n_fish_fleets` specifying the estimation
  structure for annual selectivity deviations. Required when any fleet
  has continuous time-variation. See
  [`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md)
  for all options including age-sharing options for semi-parametric
  forms.

- corr_opt_semipar:

  Character vector of length `n_fish_fleets` controlling which
  correlation components to suppress in semi-parametric (3D GMRF or 2D
  AR1) time-varying selectivity. Set to `NA` (default) for no
  suppression. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md)
  for valid suppression codes. Cohort-correlation options are invalid
  for `"2dar1"`.

- Use_fish_selex_prior:

  Integer flag. `1` = apply lognormal priors to selectivity parameters;
  `0` = no priors (default). Requires `fish_selex_prior`.

- fish_selex_prior:

  Data frame of selectivity prior hyperparameters, one row per prior.
  Required columns: `region`, `fleet`, `block`, `sex`, `par` (parameter
  index within the functional form), `mu`, `sd`, plus an optional `type`
  giving each row's target: `"par"` (the default when the column is
  absent) is a lognormal prior on one fixed selectivity parameter, with
  `mu` on the natural scale and `sd` on the log scale; `"value"` is a
  normal prior on the realized selectivity value at one bin, with both
  on the natural scale, where `par` instead names the bin (on ages or
  lengths per `fish_selex_type`) and the value is read at the first
  model year of `block`. A `"value"` row constrains the derived
  selectivity value rather than the parameters, matching the ADMB
  convention of pinning selectivity at a reference age near one, which
  no set of independent parameter priors can express. Only used when
  `Use_fish_selex_prior = 1`.

- Use_fish_selex_penalty:

  Integer (0/1). Whether a centering penalty is applied to sets of
  fishery selectivity fixed-effect parameters. Default `0`.

- fish_sel_norm_bins:

  List with one element per fishery fleet naming the bins the mean-one
  standardization averages over, or `NULL` for fleets standardizing over
  every bin (`Selex_Model = 9` only). A gear whose catchability is
  defined against part of the bin range standardizes over that part, and
  catchability absorbs the difference in scale. Default `NULL`.

- fish_sel_bin_dev_bins:

  List with one element per fishery fleet naming the bins that fleet
  overrides, or `NULL` for fleets with no overrides (e.g.
  `list(1, NULL)` frees bin 1 of fleet 1 only). An overridden bin takes
  a freely estimated annual value \\\exp(\epsilon\_{y,b})\\ in place of
  whatever the functional form produced, applied after every other
  transformation including standardization. The rest of the curve keeps
  its parametric shape. Default `NULL`.

- fishsel_pe_wt:

  Numeric vector of length `n_fish_fleets`. Per-fleet multiplier on the
  fishery selectivity process error likelihood. Default `1` for every
  fleet. `0` skips that fleet's process error likelihood altogether, so
  the deviations stay estimated but enter the objective only through the
  data and any explicit smoothness or centering penalties, which is how
  several existing assessments constrain them. Values other than 0 or 1
  make an estimated process error sigma reinterpretable, so prefer 0 or
  1 unless deliberately down-weighting. Applies only to
  `ln_fishsel_devs`; the bin-override deviations carry their own process
  error and are not affected.

- fishsel_rw_init_sigma:

  Numeric vector of length `n_fish_fleets`. Standard deviation given to
  the first year of an `"rw"` deviation series. Default `5`, which
  leaves that year effectively free. `NA` instead starts the walk at
  zero under the walk's own estimated sigma, making the first year as
  smooth as every later step. Appropriate when the base parametric curve
  already describes the first year well.

- cont_tv_fishsel_bin_devs:

  Character vector of length `n_fish_fleets` giving the process error on
  the bin-override deviations for each fleet: `"none"` (default),
  `"iid"`, or `"rw"`. A random walk carries its own estimated sigma per
  bin, with `fishsel_bin_devs_rw_init_sigma` governing its first year.

- fish_selex_penalty:

  Data frame of centering penalty specifications, required when
  `Use_fish_selex_penalty = 1`. Required columns: `region`, `fleet`,
  `block`, `sex`, `par`, and `wt`. Each row penalizes
  `wt * (log(mean(exp(pars))))^2` over the set of parameters named in
  `par`, which may be a single index or a list column of integer vectors
  naming a whole set. This pins the scalar of a non-parametric curve
  that catchability or fishing mortality would otherwise absorb, and is
  softer than fixing a bin outright. Intended for parameter sets held on
  the log scale. Default `NULL`.

- fishsel_devs_shared_bins:

  List of integer vectors grouping age or length bins that share a
  single deviation series. Only used when `fish_sel_devs_spec` contains
  one of the `"est_shared_b"` variants. Example:
  `list(1:5, 6:10, 11:30)`.

- fish_selex_type:

  Character scalar specifying whether selectivity is age- or
  length-based. Options:

  `"age"`

  :   Selectivity is defined over age bins.

  `"length"`

  :   Selectivity is defined over length bins.

  Determines the bin dimension used for all fishery selectivity
  functions, including parametric, time-varying, and non-parametric
  forms.

- use_fixed_fish_sel:

  Integer vector of length `n_fish_fleets` indicating whether fishery
  selectivity is fixed (`1`) or estimated (`0`) for each fleet.

- fish_sel_input:

  Array of fixed fishery selectivity values with dimensions:
  `[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]`.
  Required when any element of `use_fixed_fish_sel == 1`.

- fish_sel_nonpar_est_bins:

  Optional list defining bin groupings for non-parametric fishery
  selectivity. Structure is `[[fleet]][[block]]`, where each element is
  a list of integer vectors. Each vector defines a group of bins that
  share a single estimated selectivity parameter. Indices must
  correspond to the bin dimension defined by `fish_selex_type`.

- fish_sel_sex_offset:

  Character vector of length `n_fish_fleets` linking the sexes of a
  fleet's selectivity, for models with `n_sexes > 1`. Options per fleet:

  `"none"` (default)

  :   Each sex's stored parameters are its own, exactly as before this
      option existed.

  `"par"`

  :   The stored fixed-effect parameter slots of every sex beyond the
      first hold additive offsets on the first sex's stored
      (transformed-scale) parameters, so for log-scale parameters the
      sex-\\s\\ natural value is the first sex's times \\e^{\delta}\\.
      Offsets fixed at zero reproduce sex-shared parameters; estimating
      them links the sexes through the offset the way several existing
      assessments parameterize male selectivity.

  `"scale"`

  :   Each sex keeps its own parameters, and every sex beyond the first
      additionally carries a constant log-scale offset on the whole
      realized curve, `exp(ln_fishsel_sex_scale)`, estimated per region,
      block, and sex. The scaled curve may exceed one. Refused for
      non-parametric forms and semi-parametric time variation, whose
      post-hoc standardization would cancel a constant multiplier.

  `"apical"`

  :   Each sex keeps its own parameters, and for every sex beyond the
      first the double normal builds its limbs up to
      `exp(ln_*sel_sex_scale)` rather than to one. Selectivity at the
      first and last bins stays where that sex's own parameters put it,
      so the offset moves the middle of the curve and leaves its ends
      anchored. Requires the double normal.

  `"par_apical"`

  :   Both a par offset and an apical offset.

  `"par_scale"`

  :   Both a par offset and a scale offset.

- fish_sel_dbnrml_raw:

  `NULL` (default) or a 0/1 matrix `[n_fish_fleets x 2]` for fleets on
  the double normal: column one leaves the ascending limb as a raw
  Gaussian instead of anchoring it to `p5` at the first bin, column two
  does the same for the descending limb and `p6`.

- fish_sel_dbnrml_startbin:

  `NULL` (default) or an integer vector `[n_fish_fleets]`, the bin each
  fleet's double normal anchors its ascending limb at (`1` is the first
  bin). Bins below it take the squared ratio of their bin to it times
  the selectivity there, Stock Synthesis's convention when the
  compositions start above the population's first length bin.

- cont_tv_ret_sel:

  Character vector of length `n_fish_fleets` specifying continuous
  time-varying selectivity per fleet. Each element must be
  `"<type>_Fleet_<f>"`. Valid types:

  `"none"`

  :   No continuous time-variation (default).

  `"iid"`

  :   IID annual deviations on selectivity parameters.

  `"rw"`

  :   Random walk in selectivity parameters over time.

  `"3dmarg"`

  :   3D GMRF with marginal variance parameterization.

  `"3dcond"`

  :   3D GMRF with conditional variance parameterization.

  `"2dar1"`

  :   2D separable AR1 in bin and year dimensions.

  If any fleet has `cont_tv_ret_sel != "none"`, both
  `retsel_pe_pars_spec` and `ret_sel_devs_spec` must also be provided.

- ret_sel_blocks:

  Character vector defining discrete selectivity time blocks per fleet.
  Each element follows `"Block_<b>_Year_<s>-<e>_Fleet_<f>"` or
  `"Block_<b>_Year_<s>-terminal_Fleet_<f>"`. Use `"none_Fleet_<f>"`
  (default) for a single constant block. Blocks must be non-overlapping
  and together span all model years for the specified fleet. Mutually
  exclusive with `cont_tv_ret_sel != "none"` for the same fleet.

- ret_sel_model:

  Character vector specifying the selectivity functional form for each
  fleet (and optionally each time block). Each element must follow one
  of:

  - `"<model>_Fleet_<f>"`: single form for all years of fleet `f`.

  - `"<model>_Fleet_<f>_Block_<b>"`: form specific to block `b` of fleet
    `f`, as defined in `ret_sel_blocks`.

  Available models:

  `"logist1"`

  :   Logistic with \\a\_{50}\\ and slope \\k\\ (2 parameters).

  `"logist2"`

  :   Logistic with \\a\_{50}\\ and \\a\_{95}\\ (2 parameters).

  `"gamma"`

  :   Dome-shaped gamma with \\a\_{max}\\ and \\\delta\\ (2 parameters).

  `"exponential"`

  :   Exponential with a single power parameter (1 parameter).

  `"dbnrml"`

  :   Double-normal with 6 parameters.

  `"nonpar"`

  :   Non-parametric over discrete age or length bins, on the logit
      scale, then mean-standardized jointly over years and bins so the
      grand mean of the surface is one. Bins may be grouped through the
      non-parametric bin mapping. No fixed functional form is imposed.

  `"nonparlog"`

  :   Non-parametric on the log scale, standardized so each year's
      selectivity averages to one over `*_sel_norm_bins`. Only
      within-year contrasts are identified; the level is absorbed by
      catchability or fishing mortality.

  `"nonparfree"`

  :   Non-parametric on the log scale with no standardization,
      \\\exp(\theta)\\, so the values carry the height of the curve as
      well as its shape. This is the form for a stream fit age by age: a
      free catchability per age and a selectivity estimated at age are
      one quantity written two ways, so the whole age multiplier lives
      here and no catchability is set. Pin one bin, by leaving it out of
      the estimated bins, whenever the mean it multiplies is also free.

  `"asymplogist1"`

  :   Logistic selectivity with \\a\_{50}\\ and slope \\k\\ and
      asymptotic control (3 parameters).

  `"asymplogist2"`

  :   Logistic selectivity with with \\a\_{50}\\ and \\a\_{95}\\ and
      asymptotic control (3 parameters).

  `"bicubic"`

  :   Bicubic spline over a bin-node x year-node grid, specified as
      `"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"` (optionally
      with `_Block_k`, `_SelStyr_<year>`, and/or `_NSelBins_<n>`); see
      `fish_sel_model` above for the full syntax and
      [`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md)
      (`Selex_Model == 8`) for the underlying math.

  See the model equations vignette for mathematical definitions.

- retsel_pe_pars_spec:

  Character vector of length `n_fish_fleets` specifying the estimation
  structure for selectivity process error hyperparameters. Required when
  any fleet has continuous time-variation. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md)
  for all options.

- ret_fixed_sel_pars_spec:

  Character vector of length `n_fish_fleets` specifying how fixed-effect
  selectivity parameters are estimated. See
  [`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md)
  for all options (`"est_all"`, `"est_shared_r"`, `"est_shared_s"`,
  `"est_shared_r_s"`, `"est_shared_f_x"`, `"fix"`).

- ret_sel_devs_spec:

  Character vector of length `n_fish_fleets` specifying the estimation
  structure for annual selectivity deviations. Required when any fleet
  has continuous time-variation. See
  [`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md)
  for all options including age-sharing options for semi-parametric
  forms.

- ret_sel_corr_opt_semipar:

  Character vector of length `n_fish_fleets` controlling which
  correlation components to suppress in semi-parametric (3D GMRF or 2D
  AR1) time-varying selectivity. Set to `NA` (default) for no
  suppression. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md)
  for valid suppression codes. Cohort-correlation options are invalid
  for `"2dar1"`.

- Use_ret_selex_prior:

  Integer flag. `1` = apply lognormal priors to selectivity parameters;
  `0` = no priors (default). Requires `ret_selex_prior`.

- ret_selex_prior:

  Data frame of selectivity prior hyperparameters. Required columns:
  `region`, `fleet`, `block`, `sex`, `par`, `mu`, `sd`, plus an optional
  `type` (`"par"`/`"value"`; see `fish_selex_prior`).

- retsel_devs_shared_bins:

  List of integer vectors grouping age or length bins that share a
  single deviation series. Only used when `ret_sel_devs_spec` contains
  one of the `"est_shared_b"` variants. Example:
  `list(1:5, 6:10, 11:30)`.

- retsel_pe_wt:

  Numeric vector of length `n_fish_fleets`. Per-fleet multiplier on the
  retention selectivity process error likelihood, the retention
  counterpart of `fishsel_pe_wt`. Default `1` for every fleet, and `0`
  skips that fleet's process error likelihood so its deviations stay
  estimated but are constrained only by the data and any explicit
  smoothness or centering penalties.

- retsel_rw_init_sigma:

  Numeric vector of length `n_fish_fleets`. Standard deviation given to
  the first year of an `"rw"` retention deviation series, the retention
  counterpart of `fishsel_rw_init_sigma`. Default `5`; `NA` instead
  starts the walk at zero under the walk's own estimated sigma.

- ret_selex_type:

  Character scalar specifying whether retained selectivity is age- or
  length-based. Options:

  `"age"`

  :   Selectivity is defined over age bins.

  `"length"`

  :   Selectivity is defined over length bins.

  Determines the bin dimension used for all retained selectivity
  functions, including parametric, time-varying, and non-parametric
  forms.

- use_fixed_ret_sel:

  Integer vector of length `n_fish_fleets` indicating whether to fix
  selectivity (`1`) or estimate it (`0`).

- ret_sel_input:

  Array of fixed selectivity values with dimensions
  `[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]`.

- ret_sel_nonpar_est_bins:

  Optional list specifying bin groupings for non-parametric retained
  selectivity. Structure is `[[fleet]][[block]]`, where each element is
  a list of bin index vectors defining grouped parameters.

- ret_sel_sex_offset:

  Character vector of length `n_fish_fleets` linking the sexes of a
  fleet's retention curve, with the options and meaning of
  `fish_sel_sex_offset`. Default `"none"`. Retention is a fraction, so a
  scale offset is only sensible where the scaled curve stays at or below
  one.

- ...:

  Optional starting value overrides for selectivity parameters.

## Value

The input `input_list` with `$data`, `$par`, and `$map` updated. Key
additions include the parsed integer arrays for `cont_tv_fish_sel`,
`fish_sel_blocks`, `fish_sel_model`, and `fish_q_blocks`; starting value
arrays for all four parameter groups; and their corresponding factor
maps.

## Details

Selectivity time-variation and blocked selectivity are mutually
exclusive within a fleet. Specifying both for the same fleet will raise
an error.

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
