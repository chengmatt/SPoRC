# Set up survey selectivity and catchability specifications

Configures all survey selectivity and catchability components of the
estimation model: continuous and blocked time-varying selectivity,
selectivity functional forms, catchability blocks and optional
environmental covariate effects, process error and deviation mapping,
and selectivity/catchability priors. Delegates parameter mapping to four
internal helpers
([`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md),
[`do_q_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_q_mapping.md),
[`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md),
[`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md)).
Must be called after
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)
and before model compilation.

## Usage

``` r
Setup_Mod_Srvsel_and_Q(
  input_list,
  cont_tv_srv_sel = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ""),
  srv_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ""),
  srv_sel_model,
  Use_srv_q_prior = 0,
  srv_q_prior = NA,
  srv_q_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ""),
  srvsel_pe_pars_spec = NULL,
  srv_fixed_sel_pars_spec,
  srv_q_spec = NULL,
  srv_q_type = rep("est", input_list$data$n_srv_fleets),
  srv_sel_devs_spec = NULL,
  corr_opt_semipar = NULL,
  srv_q_formula = NULL,
  srv_q_cov_dat = NULL,
  Use_srv_selex_prior = 0,
  srv_selex_prior = NULL,
  Use_srv_selex_penalty = 0,
  srv_sel_norm_bins = NULL,
  srv_sel_bin_dev_bins = NULL,
  srvsel_pe_wt = rep(1, input_list$data$n_srv_fleets),
  srvsel_rw_init_sigma = rep(5, input_list$data$n_srv_fleets),
  cont_tv_srvsel_bin_devs = rep("none", input_list$data$n_srv_fleets),
  srv_selex_penalty = NULL,
  t_srv = array(1, dim = c(input_list$data$n_regions, input_list$data$n_seas,
    input_list$data$n_srv_fleets)),
  srvsel_devs_shared_bins = NULL,
  srv_selex_type = "age",
  use_fixed_srv_sel = rep(0, input_list$data$n_srv_fleets),
  srv_sel_input = NULL,
  srv_sel_nonpar_est_bins = NULL,
  srv_sel_sex_offset = rep("none", input_list$data$n_srv_fleets),
  srv_sel_dbnrml_raw = NULL,
  srv_sel_dbnrml_startbin = NULL,
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions. `$data$srv_selex_type` must
  already be set by
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

- cont_tv_srv_sel:

  Character vector defining the continuous time-variation form per
  fleet. Each element follows `"<type>_Fleet_x"`. Options:

  `"none"`

  :   No continuous time variation (default).

  `"iid"`

  :   IID deviations across years.

  `"rw"`

  :   Random walk in time.

  `"3dmarg"`

  :   3D marginal GMRF (age × year × cohort).

  `"3dcond"`

  :   3D conditional GMRF.

  `"2dar1"`

  :   2D AR1 (bin × year).

  When any fleet uses a non-`"none"` type, both `srvsel_pe_pars_spec`
  and `srv_sel_devs_spec` must be specified. Default: `"none_Fleet_x"`
  for each fleet.

- srv_sel_blocks:

  Character vector defining discrete time blocks for survey selectivity.
  Each element follows `"Block_k_Year_a-b_Fleet_x"` or `"none_Fleet_x"`
  (constant selectivity). Use `"terminal"` in place of the end year to
  extend to the final model year. Parsed into an internal array
  `[n_regions × n_years × n_srv_fleets]`. Blocked and continuous
  time-varying selectivity are mutually exclusive for a given fleet.
  Default: `"none_Fleet_x"` for each fleet.

- srv_sel_model:

  Character vector specifying the selectivity functional form per fleet,
  and optionally per time block. Each element follows one of:

  `"<model>_Fleet_x"`

  :   Single form applied across all years for fleet `x`.

  `"<model>_Fleet_x_Block_k"`

  :   Form applied only to block `k` for fleet `x`, as defined in
      `srv_sel_blocks`. Required when multiple blocks are defined for a
      fleet.

  Available models (see the model equations vignette for
  parameterisations):

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
      defined via `srv_sel_blocks`). An optional `_SelStyr_<year>`
      suffix (a calendar year within the block) restricts the actual
      spline fit to `SelStyr`:block-end; years within the block before
      `SelStyr` are held constant at the `SelStyr` year's fitted curve,
      rather than fitting the surface over the whole block. An optional
      `_NSelBins_<n>` suffix restricts the actual spline fit to the
      first `n` bins (ages or lengths, per `srv_selex_type`); bins
      beyond `n` are held constant at the last fitted bin's curve.

  No default; must be provided.

- Use_srv_q_prior:

  Integer (0/1). Whether log-normal priors are applied to survey
  catchability parameters. Default `0`.

- srv_q_prior:

  Data frame of catchability prior specifications. Required columns:
  `region`, `fleet`, `block`, `mu` (prior mean on natural scale), `sd`
  (prior SD on log scale). Each row specifies a
  \\\log\text{N}(\log(\mu), \text{sd})\\ prior. Ignored when
  `Use_srv_q_prior = 0`. Default `NA`.

- srv_q_blocks:

  Character vector defining discrete time blocks for survey
  catchability. Same format as `srv_sel_blocks`:
  `"Block_k_Year_a-b_Fleet_x"` or `"none_Fleet_x"`. Parsed into an array
  `[n_regions × n_years × n_srv_fleets]`. Default: `"none_Fleet_x"` for
  each fleet.

- srvsel_pe_pars_spec:

  Character vector `[n_srv_fleets]` or `NULL`. Sharing structure for
  process error hyperparameters. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md)
  for full option descriptions. Default `NULL`.

- srv_fixed_sel_pars_spec:

  Character vector `[n_srv_fleets]`. Sharing structure for fixed-effect
  selectivity parameters. See
  [`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md)
  for full option descriptions. No default; must be provided.

- srv_q_spec:

  Character vector `[n_srv_fleets]` or `NULL`. Sharing structure for
  catchability. See
  [`do_q_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_q_mapping.md)
  for full option descriptions. Default `NULL`.

- srv_q_type:

  Character vector `[n_srv_fleets]` controlling how catchability is
  obtained. `"est"` (default) estimates `ln_srv_q`. `"arith"`
  concentrates it out of the likelihood as the ratio of mean observed to
  mean predicted index, and `"geo"` does the same on the log scale as
  `exp(mean(log(obs) - log(pred)))`. Both analytic forms solve one
  catchability per region and fleet using only the years with
  observations, ignore any block structure, and fix that fleet's
  `ln_srv_q` regardless of `srv_q_spec`.

- srv_sel_devs_spec:

  Character vector `[n_srv_fleets]` or `NULL`. Sharing structure for
  selectivity deviation time series. See
  [`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md)
  for full option descriptions. Default `NULL`.

- corr_opt_semipar:

  Character vector `[n_srv_fleets]` or `NULL`. Specifies correlation
  components to suppress for 3D GMRF or 2D AR1 forms. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md)
  for valid values. Default `NULL`.

- srv_q_formula:

  Named list of R formulas specifying environmental covariate
  relationships for catchability per region-fleet combination. Names
  follow the convention `"Region_r_Fleet_f"`. Covariates must be present
  in `srv_q_cov_dat`. If `NULL`, no covariate effects are included.
  Default `NULL`.

- srv_q_cov_dat:

  Named list of numeric vectors (length = `n_years`) containing
  covariate time series referenced in `srv_q_formula`. All vectors must
  be the same length and contain no missing values; set values to `0`
  for years when the survey is not active. If `NULL`, covariate effects
  are excluded. Default `NULL`.

- Use_srv_selex_prior:

  Integer (0/1). Whether log-normal priors are applied to survey
  selectivity parameters. Default `0`.

- srv_selex_prior:

  Data frame of selectivity prior specifications, one row per prior.
  Required columns: `region`, `fleet`, `block`, `sex`, `par`, `mu`,
  `sd`, plus an optional `type` giving each row's target: `"par"` (the
  default when the column is absent) is a lognormal prior on one fixed
  selectivity parameter, with `mu` on the natural scale and `sd` on the
  log scale; `"value"` is a normal prior on the realized selectivity
  value at one bin, with both on the natural scale, where `par` instead
  names the bin (on ages or lengths per `srv_selex_type`) and the value
  is read at the first model year of `block`. A `"value"` row constrains
  the derived selectivity value rather than the parameters, matching the
  ADMB convention of pinning survey selectivity at a reference age near
  one, which no set of independent parameter priors can express. Ignored
  when `Use_srv_selex_prior = 0`. Default `NULL`.

- Use_srv_selex_penalty:

  Integer (0/1). Whether a centering penalty is applied to sets of
  survey selectivity fixed-effect parameters. Default `0`.

- srv_sel_norm_bins:

  List with one element per survey fleet naming the bins the mean-one
  standardization averages over, or `NULL` for fleets standardizing over
  every bin (`Selex_Model = 9` only). A gear whose catchability is
  defined against part of the bin range standardizes over that part, and
  catchability absorbs the difference in scale. Default `NULL`.

- srv_sel_bin_dev_bins:

  List with one element per survey fleet naming the bins that fleet
  overrides, or `NULL` for fleets with no overrides (e.g.
  `list(1, NULL)` frees bin 1 of fleet 1 only). An overridden bin takes
  a freely estimated annual value \\\exp(\epsilon\_{y,b})\\ in place of
  whatever the functional form produced, applied after every other
  transformation including standardization. The rest of the curve keeps
  its parametric shape. Default `NULL`.

- srvsel_pe_wt:

  Numeric vector `[n_srv_fleets]`. Per-fleet multiplier on the survey
  selectivity process error likelihood. Default `1` for every fleet. `0`
  skips that fleet's process error likelihood altogether, so the
  deviations stay estimated but enter the objective only through the
  data and any explicit smoothness or centering penalties, which is how
  several existing assessments constrain them. Values other than 0 or 1
  make an estimated process error sigma reinterpretable, so prefer 0 or
  1 unless deliberately down-weighting. Applies only to
  `ln_srvsel_devs`; the bin-override deviations carry their own process
  error and are not affected.

- srvsel_rw_init_sigma:

  Numeric vector `[n_srv_fleets]`. Standard deviation given to the first
  year of an `"rw"` deviation series. Default `5`, which leaves that
  year effectively free. `NA` instead starts the walk at zero under the
  walk's own estimated sigma, making the first year as smooth as every
  later step. Appropriate when the base parametric curve already
  describes the first year well.

- cont_tv_srvsel_bin_devs:

  Character vector `[n_srv_fleets]` giving the process error on the
  bin-override deviations for each fleet: `"none"` (default), `"iid"`,
  or `"rw"`. A random walk carries its own estimated sigma per bin, with
  `srvsel_bin_devs_rw_init_sigma` governing its first year.

- srv_selex_penalty:

  Data frame of centering penalty specifications, required when
  `Use_srv_selex_penalty = 1`. Required columns: `region`, `fleet`,
  `block`, `sex`, `par`, and `wt`. Each row penalizes
  `wt * (log(mean(exp(pars))))^2` over the set of parameters named in
  `par`, which may be a single index or a list column of integer vectors
  naming a whole set. This pins the scalar of a non-parametric curve
  that catchability or fishing mortality would otherwise absorb, and is
  softer than fixing a bin outright. Intended for parameter sets held on
  the log scale. Default `NULL`.

- t_srv:

  Survey timing fraction within a given year (annual models) or season
  (seasonal models), array `[n_regions × n_seas × n_srv_fleets]`.
  Default: `1` (end of period).

- srvsel_devs_shared_bins:

  List of integer vectors defining bin groups for age/length-sharing of
  deviations under semi-parametric forms (e.g.,
  `list(1:5, 6:10, 11:30)`). Required when `srv_sel_devs_spec` includes
  any `"est_shared_b"` variant. Default `NULL`.

- srv_selex_type:

  Character. Whether survey selectivity type is 'age' or 'length' based.
  Default: `age`.

- use_fixed_srv_sel:

  Integer vector of length `n_srv_fleets` indicating whether survey
  selectivity is fixed (`1`) or estimated (`0`) for each survey index.

- srv_sel_input:

  Array of fixed survey selectivity values used when
  `use_fixed_srv_sel == 1`. Dimensions:
  `[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_srv_fleets]`.
  Required whenever any survey has fixed selectivity specified.

- srv_sel_nonpar_est_bins:

  Optional list defining bin groupings for non-parametric survey
  selectivity. Structure is `[[survey]][[block]]`, where each element is
  a list of integer vectors. Each vector defines a group of bins that
  share a single estimated selectivity parameter. Indices must
  correspond to the bin dimension defined by the survey selectivity type
  (age or length).

- srv_sel_sex_offset:

  Character vector of length `n_srv_fleets` linking the sexes of a
  fleet's selectivity, for models with `n_sexes > 1`. Options per fleet
  are `"none"` (default, each sex's stored parameters are its own),
  `"par"` (sexes beyond the first store additive offsets on the first
  sex's transformed-scale parameters), `"scale"` (sexes beyond the first
  carry an estimated constant log-scale offset on the whole realized
  curve), `"par_scale"` (both), `"apical"` (sexes beyond the first build
  their double normal's limbs up to an estimated height rather than to
  one, which moves the middle of the curve and leaves its two ends where
  their own parameters put them), and `"par_apical"` (both). See
  `fish_sel_sex_offset` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)
  for the full description.

- srv_sel_dbnrml_raw:

  `NULL` (default) or a 0/1 matrix `[n_srv_fleets x 2]` for fleets on
  the double normal: column one leaves the ascending limb as a raw
  Gaussian instead of anchoring it to `p5` at the first bin, column two
  does the same for the descending limb and `p6`.

- srv_sel_dbnrml_startbin:

  `NULL` (default) or an integer vector `[n_srv_fleets]`, the bin each
  survey's double normal anchors its ascending limb at; see
  `fish_sel_dbnrml_startbin` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md).

- ...:

  Optional named starting values for selectivity and catchability
  parameters.

## Value

The input `input_list` with selectivity and catchability configuration
stored in `$data` (`cont_tv_srv_sel`, `srv_sel_blocks`, `srv_sel_model`,
`srv_q_blocks`, `srv_q_prior`, `Use_srv_q_prior`, `do_srv_q_cov`,
`srv_q_cov`, `Use_srv_selex_prior`, `srv_selex_prior`, `t_srv`);
starting values in `$par` for `srv_fixed_sel_pars`, `ln_srv_q`,
`srvsel_pe_pars`, `ln_srvsel_devs`, and `srv_q_coeff`; and factor maps
in `$map` for all five parameter arrays plus `srv_q_coeff`.

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
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
