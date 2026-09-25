# Set up retained fishery selectivity

Sets the retention selectivity forms, time blocks, continuous time
variation, process error, annual deviations and fixed or estimated
parameters. Time variation through `cont_tv_ret_sel` and blocks through
`ret_sel_blocks` are mutually exclusive within a fleet. Call after
[`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md).

## Usage

``` r
Setup_Mod_Retsel(
  input_list,
  cont_tv_ret_sel,
  ret_sel_blocks,
  ret_sel_model,
  retsel_pe_pars_spec,
  ret_fixed_sel_pars_spec,
  ret_sel_devs_spec,
  ret_sel_corr_opt_semipar,
  Use_ret_selex_prior,
  ret_selex_prior,
  retsel_devs_shared_bins,
  ret_selex_type,
  use_fixed_ret_sel,
  ret_sel_input,
  ret_sel_bin_dev_bins = NULL,
  cont_tv_retsel_bin_devs = rep("none", input_list$data$n_fish_fleets),
  retsel_pe_wt = rep(1, input_list$data$n_fish_fleets),
  retsel_rw_init_sigma = rep(5, input_list$data$n_fish_fleets),
  retsel_dont_est_dev_first = rep(0, input_list$data$n_fish_fleets),
  ret_sel_nonpar_est_bins,
  ret_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- cont_tv_ret_sel:

  Character vector of continuous time variation per fleet, each
  `"<type>_Fleet_<f>"`: `"none"`, `"iid"`, `"rw"`, `"3dmarg"`,
  `"3dcond"` or `"2dar1"`. Stored as integer codes.

- ret_sel_blocks:

  Character vector of discrete blocks, either `"none_Fleet_<f>"` or
  `"Block_<b>_Year_<start>-<end>_Fleet_<f>"` with `"terminal"` allowed
  as the end year. Expanded into an
  `[n_regions × n_years × n_fish_fleets]` array.

- ret_sel_model:

  Character vector of the retention form per fleet, and optionally per
  block: `"<type>_Fleet_<f>"` or `"<type>_Fleet_<f>_Block_<b>"`. The
  forms and their syntax are those of `fish_sel_model` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md).

- retsel_pe_pars_spec:

  Process error parameters for time-varying retention, one entry per
  fishery fleet.

- ret_fixed_sel_pars_spec:

  Which fixed retention parameters are estimated.

- ret_sel_devs_spec:

  Structure of the annual retention deviations per fleet.

- ret_sel_corr_opt_semipar:

  Optional correlation structure for semi-parametric retention
  deviations, one entry per fishery fleet.

- Use_ret_selex_prior:

  Integer flag (0/1) for retention selectivity priors.

- ret_selex_prior:

  Data frame with columns `region`, `fleet`, `block`, `sex`, `par`,
  `mu`, `sd` and an optional `type` (`"par"` or `"value"`), as
  `fish_selex_prior` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md).

- retsel_devs_shared_bins:

  Bins sharing one retention deviation series.

- ret_selex_type:

  Character scalar, `"age"` or `"length"`.

- use_fixed_ret_sel:

  Binary array `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`,
  `1` to use `ret_sel_input`.

- ret_sel_input:

  Fixed retention array
  `[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]`.

- ret_sel_bin_dev_bins:

  List with one element per fleet naming the bins that fleet overrides,
  or `NULL` for none. An overridden bin takes a freely estimated annual
  value in place of what the functional form produced, applied after
  every other transformation including standardization. Default `NULL`.

- cont_tv_retsel_bin_devs:

  Character vector `[n_fish_fleets]` of the process error on the
  bin-override deviations: `"none"` (default), `"iid"` or `"rw"`.

- retsel_pe_wt:

  Numeric vector `[n_fish_fleets]` multiplying the retention process
  error likelihood. Default `1`. `0` skips that fleet's process error,
  so the deviations stay estimated but enter the objective only through
  the data and any smoothness or centering penalties. Values other than
  0 or 1 make an estimated process error sigma reinterpretable. Applies
  to `ln_retsel_devs` only.

- retsel_rw_init_sigma:

  Numeric vector `[n_fish_fleets]` giving the standard deviation of the
  first year of an `"rw"` deviation series. Default `5`, which leaves
  that year effectively free. `NA` instead starts the walk at zero under
  the walk's own estimated sigma.

- ret_sel_nonpar_est_bins:

  Estimated bins for non-parametric retention.

- ret_sel_sex_offset:

  Character vector `[n_fish_fleets]` linking the sexes of a fleet's
  retention curve: `"none"` (default), `"par"`, `"scale"` or
  `"par_scale"`, as `fish_sel_sex_offset` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md).
  Retention is a fraction, so a scale offset only makes sense where the
  scaled curve stays at or below one, and nothing enforces that.

- ...:

  Optional starting values for the selectivity parameters and
  deviations.

## Value

`input_list` with the parsed retention structure arrays in `$data`, the
starting values in `$par` and the factor maps in `$map`.

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
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
