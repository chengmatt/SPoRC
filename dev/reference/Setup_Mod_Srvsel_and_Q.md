# Set up survey selectivity and catchability specifications

Sets the survey selectivity forms, time blocks, continuous time
variation, process error and deviations, the catchability blocks and
estimation structure, and the selectivity and catchability priors. Call
after
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md).

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
  srv_q_model = NULL,
  sigma_srv_q_spec = "est_all",
  srv_q_rho_spec = "est_all",
  srv_q_rw_init_sigma = NA,
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
  srvsel_dont_est_dev_first = rep(0, input_list$data$n_srv_fleets),
  srv_sel_sex_offset = rep("none", input_list$data$n_srv_fleets),
  srv_sel_dbnrml_raw = NULL,
  srv_sel_dbnrml_startbin = NULL,
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`.
  `$data$srv_selex_type` must already be set by
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

- cont_tv_srv_sel:

  Character vector of continuous time variation per fleet, each
  `"<type>_Fleet_x"`: `"none"` (default), `"iid"`, `"rw"`, `"3dmarg"`,
  `"3dcond"` (3D GMRF over age, year and cohort) or `"2dar1"` (separable
  over bin and year). Any fleet other than `"none"` also needs
  `srvsel_pe_pars_spec` and `srv_sel_devs_spec`.

- srv_sel_blocks:

  Character vector of discrete selectivity time blocks, each
  `"Block_k_Year_a-b_Fleet_x"` with `"terminal"` allowed as the end
  year, or `"none_Fleet_x"` (default) for constant selectivity. Parsed
  into an `[n_regions × n_years × n_srv_fleets]` array. Mutually
  exclusive with continuous time variation for the same fleet.

- srv_sel_model:

  Character vector of the selectivity form per fleet, and optionally per
  block: `"<model>_Fleet_x"` or `"<model>_Fleet_x_Block_k"`, the latter
  required when a fleet has several blocks. The forms and their syntax,
  including the `"bicubic"` suffixes, are those of `fish_sel_model` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md).
  No default.

- Use_srv_q_prior:

  Integer (0/1) for lognormal priors on survey catchability. Default
  `0`.

- srv_q_prior:

  Data frame with columns `region`, `fleet`, `block`, `mu` on the
  natural scale and `sd` on the log scale, one row per
  \\\log\text{N}(\log(\mu), \text{sd})\\ prior. Default `NA`.

- srv_q_blocks:

  Character vector of discrete catchability time blocks, in the format
  of `srv_sel_blocks`. Parsed into an
  `[n_regions × n_years × n_srv_fleets]` array. Default
  `"none_Fleet_x"`.

- srvsel_pe_pars_spec:

  Character vector `[n_srv_fleets]` or `NULL` (default) of the sharing
  structure for the process error hyperparameters. See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md).

- srv_fixed_sel_pars_spec:

  Character vector `[n_srv_fleets]` of the sharing structure for the
  fixed-effect selectivity parameters. No default. See
  [`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md).

- srv_q_spec:

  Character vector `[n_srv_fleets]` or `NULL` (default) of the sharing
  structure for catchability. See
  [`do_q_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_q_mapping.md).

- srv_q_type:

  Character vector `[n_srv_fleets]` of how catchability is obtained.
  `"est"` (default) estimates `ln_srv_q`, `"arith"` concentrates it out
  as the ratio of mean observed to mean predicted index, and `"geo"`
  does the same on the log scale as `exp(mean(log(obs) - log(pred)))`.
  Both analytic forms use the years with observations only and fix that
  fleet's `ln_srv_q` whatever `srv_q_spec` says. The solve runs within
  each `srv_q_blocks` block.

- srv_sel_devs_spec:

  Character vector `[n_srv_fleets]` or `NULL` (default) of the sharing
  structure for the deviation series. See
  [`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md).

- corr_opt_semipar:

  Character vector `[n_srv_fleets]` or `NULL` (default) of which
  correlation components to suppress under the 3D GMRF or 2D AR1 forms.
  See
  [`do_sel_pe_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_pe_pars_mapping.md).

- srv_q_model:

  Character vector `[n_srv_fleets]` of the process error on annual
  catchability deviations: `"none"` (default), `"iid"`, `"rw"`, `"ar1"`
  or `"dsem"`, which hands the series to
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md).
  Catchability is then \\\exp(\ln q\_{r,b,f} + \epsilon\_{r,y,f})\\. A
  fleet with deviations cannot also have `srv_q_blocks` or an
  analytically solved `srv_q_type`.

- sigma_srv_q_spec:

  Sharing string for the deviation standard deviation over region and
  fleet: `"est_all"` (default), `"est_shared_r"`, `"est_shared_f"`,
  `"est_shared_r_f"` or `"fix"`.

- srv_q_rho_spec:

  Sharing string for the AR1 correlation, with the same options as
  `sigma_srv_q_spec`. Default `"est_all"`. Only read under
  `srv_q_model = "ar1"`.

- srv_q_rw_init_sigma:

  Standard deviation of the first estimated year of a random walk. `NA`
  (default) starts the walk at zero under its own sigma, which keeps
  `ln_srv_q` as the level of the series; a wide value leaves the level
  free and confounds it with `ln_srv_q`.

- Use_srv_selex_prior:

  Integer (0/1) for priors on the survey selectivity parameters. Default
  `0`.

- srv_selex_prior:

  Data frame with columns `region`, `fleet`, `block`, `sex`, `par`,
  `mu`, `sd` and an optional `type`. `"par"` (the default) is a
  lognormal prior on one fixed selectivity parameter, with `mu` on the
  natural scale and `sd` on the log scale. `"value"` is a normal prior
  on the realized selectivity at one bin, both on the natural scale,
  where `par` names the bin and the value is read at the first model
  year of `block`; that is the ADMB convention of pinning survey
  selectivity at a reference age near one, which no set of independent
  parameter priors can express. Default `NULL`.

- Use_srv_selex_penalty:

  Integer (0/1). Whether a centering penalty is applied to sets of
  survey selectivity fixed-effect parameters. Default `0`.

- srv_sel_norm_bins:

  List with one element per fleet naming the bins the mean-one
  standardization averages over, or `NULL` for fleets standardizing over
  every bin. Read under `"nonparlog"` only. A gear whose catchability is
  defined against part of the bin range standardizes over that part, and
  catchability absorbs the difference in scale. Default `NULL`.

- srv_sel_bin_dev_bins:

  List with one element per fleet naming the bins that fleet overrides,
  or `NULL` for none, e.g. `list(1, NULL)`. An overridden bin takes a
  freely estimated annual value \\\exp(\epsilon\_{y,b})\\ in place of
  what the functional form produced, applied after every other
  transformation including standardization, while the rest of the curve
  keeps its parametric shape. Default `NULL`.

- srvsel_pe_wt:

  Numeric vector `[n_srv_fleets]` multiplying the survey selectivity
  process error likelihood. Default `1`. `0` skips that fleet's process
  error, so the deviations stay estimated but enter the objective only
  through the data and any smoothness or centering penalties. Values
  other than 0 or 1 make an estimated process error sigma
  reinterpretable. Applies to `ln_srvsel_devs` only; the bin-override
  deviations have their own process error.

- srvsel_rw_init_sigma:

  Numeric vector `[n_srv_fleets]` giving the standard deviation of the
  first year of an `"rw"` deviation series. Default `5`, which leaves
  that year effectively free. `NA` instead starts the walk at zero under
  the walk's own estimated sigma, which suits a base curve that already
  describes the first year well.

- cont_tv_srvsel_bin_devs:

  Character vector `[n_srv_fleets]` of the process error on the
  bin-override deviations: `"none"` (default), `"iid"` or `"rw"`. A walk
  has its own estimated sigma per bin.

- srv_selex_penalty:

  Data frame of centering penalties with columns `region`, `fleet`,
  `block`, `sex`, `par` and `wt`, required when
  `Use_srv_selex_penalty = 1`. Each row penalizes
  `wt * (log(mean(exp(pars))))^2` over the parameters named in `par`, a
  single index or a list column of integer vectors. This pins the scalar
  of a non-parametric curve that catchability would otherwise absorb,
  and is softer than fixing a bin. Meant for parameter sets on the log
  scale. Default `NULL`.

- t_srv:

  Survey timing as a fraction of the year (annual models) or the season
  (seasonal models), array `[n_regions × n_seas × n_srv_fleets]`.
  Default `1`, the end of the period.

- srvsel_devs_shared_bins:

  List of integer vectors grouping the bins that share one deviation
  series, e.g. `list(1:5, 6:10, 11:30)`. Required when
  `srv_sel_devs_spec` names an `"est_shared_b"` variant. Default `NULL`.

- srv_selex_type:

  Character scalar, `"age"` (default) or `"length"`.

- use_fixed_srv_sel:

  Integer vector `[n_srv_fleets]`, `1` to fix survey selectivity and `0`
  to estimate it.

- srv_sel_input:

  Array of fixed survey selectivity values
  `[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_srv_fleets]`,
  required whenever any survey has fixed selectivity.

- srv_sel_nonpar_est_bins:

  Optional bin groupings for non-parametric survey selectivity,
  structured `[[survey]][[block]]`, each element a list of integer
  vectors naming the bins that share one estimated parameter. Indices
  are on the bin dim the survey selectivity type names.

- srvsel_dont_est_dev_first:

  Integer vector `[n_srv_fleets]` of 0/1, default `0`. Where `1`, that
  fleet's deviations start in year two and the fixed parameters hold
  year one. A non-parametric form has one free base parameter per bin,
  so year one's deviation is that same value written twice with only
  `srvsel_rw_init_sigma` between them, a prior on a level usually meant
  to be free. Refused for the GMRF and 2D AR1 forms, whose deviations
  are a field over years and bins rather than a walk anchored at year
  one.

- srv_sel_sex_offset:

  Character vector `[n_srv_fleets]` linking the sexes of a fleet's
  selectivity when `n_sexes > 1`: `"none"` (default), `"par"`,
  `"scale"`, `"par_scale"`, `"apical"` or `"par_apical"`. See
  `fish_sel_sex_offset` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)
  for what each one does.

- srv_sel_dbnrml_raw:

  `NULL` (default) or a 0/1 matrix `[n_srv_fleets x 2]` for fleets on
  the double normal: column one leaves the ascending limb a raw Gaussian
  instead of anchoring it to `p5` at the first bin, column two does the
  same for the descending limb and `p6`.

- srv_sel_dbnrml_startbin:

  `NULL` (default) or an integer vector `[n_srv_fleets]`, the bin each
  survey's double normal anchors its ascending limb at. See
  `fish_sel_dbnrml_startbin` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md).

- ...:

  Optional starting values for the selectivity and catchability
  parameters.

## Value

`input_list` with the selectivity and catchability configuration in
`$data` (`cont_tv_srv_sel`, `srv_sel_blocks`, `srv_sel_model`,
`srv_q_blocks`, `srv_q_prior`, `Use_srv_q_prior`, `srv_q_model`,
`Use_srv_selex_prior`, `srv_selex_prior`, `t_srv`), the starting values
in `$par` for `srv_fixed_sel_pars`, `ln_srv_q`, `srvsel_pe_pars`,
`ln_srvsel_devs`, `ln_srv_q_devs`, `ln_sigma_srv_q` and `srv_q_rho`, and
their factor maps in `$map`.

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
