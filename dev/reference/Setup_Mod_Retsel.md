# Set up retained fishery selectivity

Configures all aspects of retained fishery selectivity and catchability
for the estimation model, including functional forms, time blocks,
continuous time-varying selectivity, process error structure, annual
deviations, and fixed or estimated selectivity parameters. Must be
called after
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
  ret_sel_nonpar_est_bins,
  ret_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
  ...
)
```

## Arguments

- input_list:

  Named list containing `$data`, `$par`, `$map`, and `$verbose` sublists
  as created by upstream setup functions.

- cont_tv_ret_sel:

  Character vector defining continuous time-varying selectivity
  structure per fleet. Format: `"<type>_Fleet_<f>"` where `type` is one
  of:

  `"none"`

  :   No time variation

  `"iid"`

  :   Independent year effects

  `"rw"`

  :   Random walk process

  `"3dmarg"`

  :   3D marginal structure

  `"3dcond"`

  :   3D conditional structure

  `"2dar1"`

  :   2D AR1 structure

  Output is mapped to integer codes internally.

- ret_sel_blocks:

  Character vector defining discrete selectivity blocks. Format:
  `"none_Fleet_<f>"` or `"Block_<b>_Year_<start>-<end>_Fleet_<f>"` (or
  terminal year version). Values are expanded into a 3D array:
  `[n_regions × n_years × n_fish_fleets]`.

- ret_sel_model:

  Character vector defining selectivity functional forms. Format:
  `"<type>_Fleet_<f>"` or `"<type>_Fleet_<f>_Block_<b>"`. Supported
  types:

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

  :   Bicubic spline over a bin-node x year-node grid; see
      `ret_sel_model` in
      [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)
      for full syntax.

- retsel_pe_pars_spec:

  Specification of process error parameters for time-varying
  selectivity. Length must equal `n_fish_fleets`.

- ret_fixed_sel_pars_spec:

  Specification controlling which fixed selectivity parameters are
  estimated vs fixed.

- ret_sel_devs_spec:

  Specification of annual selectivity deviations structure per fleet.

- ret_sel_corr_opt_semipar:

  Optional correlation structure for semi-parametric selectivity
  deviations. Length must equal `n_fish_fleets`.

- Use_ret_selex_prior:

  Integer flag (`0/1`) indicating whether selectivity priors are used.

- ret_selex_prior:

  Data frame of priors for selectivity parameters. Must include columns:
  `region`, `fleet`, `block`, `sex`, `par`, `mu`, `sd`, plus an optional
  `type` (`"par"`/`"value"`; see `fish_selex_prior` in
  `Setup_Mod_Fishsel_and_Q`).

- retsel_devs_shared_bins:

  Vector defining shared bins for selectivity deviation estimation
  (e.g., age or length grouping structure).

- ret_selex_type:

  Character string indicating selectivity domain: `"age"` or `"length"`.

- use_fixed_ret_sel:

  Binary indicator array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`. `1` = use
  fixed selectivity input; `0` = estimate.

- ret_sel_input:

  Fixed selectivity input array:
  `[n_pop × n_regions × n_years × n_seas × (ages or lengths) × n_sexes × n_fish_fleets]`.

- ret_sel_bin_dev_bins:

  List with one element per fishery fleet naming the bins that fleet
  overrides, or `NULL` for fleets with no overrides. An overridden bin
  takes a freely estimated annual value in place of whatever the
  functional form produced, applied after every other transformation
  including standardization. Default `NULL`.

- cont_tv_retsel_bin_devs:

  Character vector of length `n_fish_fleets` giving the process error on
  the bin-override deviations for each fleet: `"none"` (default),
  `"iid"`, or `"rw"`.

- retsel_pe_wt:

  Numeric vector of length `n_fish_fleets`. Per-fleet multiplier on the
  retention selectivity process error likelihood. Default `1` for every
  fleet. `0` skips that fleet's process error likelihood altogether, so
  the deviations stay estimated but enter the objective only through the
  data and any explicit smoothness or centering penalties. Values other
  than 0 or 1 make an estimated process error sigma reinterpretable, so
  prefer 0 or 1 unless deliberately down-weighting. Applies only to
  `ln_retsel_devs`; the bin-override deviations carry their own process
  error and are not affected.

- retsel_rw_init_sigma:

  Numeric vector of length `n_fish_fleets`. Standard deviation given to
  the first year of an `"rw"` deviation series. Default `5`, which
  leaves that year effectively free. `NA` instead starts the walk at
  zero under the walk's own estimated sigma, making the first year as
  smooth as every later step. Appropriate when the base parametric curve
  already describes the first year well.

- ret_sel_nonpar_est_bins:

  Vector defining estimated bins for non-parametric selectivity.

- ret_sel_sex_offset:

  Character vector of length `n_fish_fleets` linking the sexes of a
  fleet's retention curve: `"none"` (default), `"par"`, `"scale"`, or
  `"par_scale"`, with the same meaning as `fish_sel_sex_offset` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md).
  Retention is a fraction, so a scale offset is only sensible where the
  scaled curve stays at or below one (a negative offset, the
  less-retained sex); nothing enforces that.

- ...:

  Optional starting values for selectivity parameters and deviations.

## Value

The updated `input_list` with:

- Parsed selectivity structure arrays in `$data`

- Starting parameter arrays in `$par`

- Factor mapping objects in `$map`

## Details

Selectivity time-variation via `cont_tv_ret_sel` and blocked selectivity
via `ret_sel_blocks` are mutually exclusive within a fleet. Specifying
both for the same fleet will result in an error.

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
