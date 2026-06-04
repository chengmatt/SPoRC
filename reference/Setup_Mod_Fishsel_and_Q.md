# Setup fishery selectivity and catchability specifications

Setup fishery selectivity and catchability specifications

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
  fishsel_pe_pars_spec = NULL,
  fish_fixed_sel_pars_spec = NULL,
  fish_q_spec = NULL,
  fish_sel_devs_spec = NULL,
  corr_opt_semipar = NULL,
  Use_fish_selex_prior = 0,
  fish_selex_prior = NULL,
  cont_tv_fish_sel_penalty = TRUE,
  fishsel_devs_shared_ages = NULL,
  ...
)
```

## Arguments

- input_list:

  List containing a data list, parameter list, and map list

- cont_tv_fish_sel:

  Character vector specifying the form of continuous time-varying
  selectivity for each fishery fleet. The vector must be length
  `n_fish_fleets`, and each element must follow the structure:
  `"<time variation type>_Fleet_<fleet number>"`.

  Valid time variation types include:

  - `"none"`: No continuous time variation (default)

  - `"iid"`: Independent and identically distributed deviations across
    years.

  - `"rw"`: Random walk in time.

  - `"3dmarg"`: 3D marginal time-varying selectivity.

  - `"3dcond"`: 3D conditional time-varying selectivity.

  - `"2dar1"`: Two-dimensional AR1 process.

  For example:

  - `"iid_Fleet_1"` applies an iid time-varying structure to Fleet 1.

  - `"none_Fleet_2"` means no time variation is used for Fleet 2.

- fish_sel_blocks:

  Character vector specifying the fishery selectivity blocks for each
  region and fleet.

  Each element must follow one of the following structures:

  - \`"Block\_\<block number\>\_Year\_\<start\>-\<end\>\_Fleet\_\<fleet
    number\>"\`

  - \`"Block\_\<block number\>\_Year\_\<start\>-terminal_Fleet\_\<fleet
    number\>"\`

  - \`"none_Fleet\_\<fleet number\>"\`

  This argument defines how fishery selectivity varies over time for
  each fleet:

  - `"Block_..."` entries specify discrete time blocks during which
    selectivity parameters are assumed constant.

  - `"none_..."` entries indicate that selectivity is constant across
    all years for the specified fleet.

  If time-block-based selectivity is specified for a fleet (via
  `fish_sel_blocks`), its corresponding continuous selectivity option
  (in `cont_tv_fish_sel`) must be set to `"none_Fleet_<fleet number>"`.
  The two approaches—blocked and continuous time-varying selectivity—are
  mutually exclusive. The default for each fleet is `"none_Fleet_x"`
  (i.e., no selectivity blocks).

- fish_sel_model:

  Character vector specifying the fishery selectivity functional form
  for each fleet, and optionally by time block.

  Each element must follow one of the following structures:

  - `"<selectivity model>_Fleet_<fleet number>"`

  - `"<selectivity model>_Fleet_<fleet number>_Block_<block number>"`

  The first form applies a single selectivity model across all years for
  the specified fleet. The second form allows the user to assign a
  distinct selectivity model to a specific time block, as defined in
  `fish_sel_blocks`.

  Available selectivity model types include:

  - `"logist1"` — Logistic function with parameters `a50` and `k`.

  - `"logist2"` — Logistic function with parameters `a50` and `a95`.

  - `"gamma"` — Dome-shaped gamma function with parameters `amax` and
    `delta`.

  - `"exponential"` — Exponential function with a power parameter.

  - `"dbnrml"` — Double-normal function with six parameters.

  If multiple selectivity time blocks are specified for a fleet (using
  `fish_sel_blocks`), then the corresponding selectivity model for each
  block must be explicitly defined using the
  `"<model>_Block_<block>_Fleet_<fleet>"` format. If blocks are not
  defined for a fleet, use the `"<model>_Fleet_<fleet number>"` format
  only. For mathematical definitions and implementation details of each
  selectivity form, refer to the model equations vignette.

- Use_fish_q_prior:

  Integer (0 or 1). Flag to enable/disable fishery catchability priors.
  When set to 1, applies log-normal priors to fishery selectivity
  parameters as specified in `fish_q_prior`. When set to 0, no priors
  are applied.

- fish_q_prior:

  Data frame containing prior specifications for fishery catchability
  parameters. Must include columns: `region` (region index), `fleet`
  (fleet index), `block` (time block index), `mu` (prior mean on natural
  scale), and `sd` (prior standard deviation on log scale). Each row
  specifies a log-normal prior N(log(mu), sd) for a given catchability
  parameter. Only parameters with rows in this data frame will have
  priors applied.

- fish_q_blocks:

  Character vector specifying fishery catchability (q) blocks for each
  fleet. Each element must follow the structure:
  `"Block_<block number>_Year_<start>-<end>_Fleet_<fleet number>"` or
  `"none_Fleet_<fleet number>"`. Default is "none_Fleet_x".

  This allows users to define time-varying catchability blocks
  independently of selectivity blocks. The blocks must be
  non-overlapping and sequential in time within each fleet.

  For example:

  - `"Block_1_Year_1-35_Fleet_1"` assigns block 1 to Fleet 1 for years
    1–35.

  - `"Block_2_Year_36-56_Fleet_1"` continues with block 2 for years
    36–56.

  - `"Block_3_Year_57-terminal_Fleet_1"` assigns block 3 from year 57 to
    the terminal year for Fleet 1.

  - `"none_Fleet_2"` indicates no catchability blocks are used for Fleet
    2.

  Internally, these specifications are converted to a
  `[n_regions, n_years, n_fish_fleets]` array, where each block is
  mapped to the appropriate years and fleets.

- fishsel_pe_pars_spec:

  Character string specifying how process error parameters for fishery
  selectivity are estimated across regions and sexes. This is only
  relevant if `cont_tv_fish_sel` is not set to `"none"`; otherwise, all
  process error parameters are treated as fixed.

  Available options include:

  - `"est_all"`: Estimates separate process error parameters for each
    region and sex.

  - `"est_shared_r"`: Shares process error parameters across regions
    (sex-specific parameters are still estimated).

  - `"est_shared_s"`: Shares process error parameters across sexes
    (region-specific parameters are still estimated).

  - `"est_shared_r_s"`: Shares process error parameters across both
    regions and sexes, estimating a single set of parameters.

  - `"est_shared_f_x"`: Shares process error parameters with another
    fleet, where `x` is the fleet number to share with. This option
    forces multiple fleets to have identical process error variance and
    correlation structures for their time-varying selectivity. For
    example, `"est_shared_f_2"` means the current fleet will use the
    same process error parameters as fleet 2. The reference fleet
    (fleet x) must use one of the other sharing options and cannot
    itself be sharing with another fleet.

  - `"fix"` or `"none"`: Does not estimate process error parameters; all
    are treated as fixed.

- fish_fixed_sel_pars_spec:

  Character string specifying the structure for estimating fixed-effect
  parameters of the fishery selectivity model (e.g., a50, k, amax). This
  controls whether selectivity parameters are estimated separately or
  shared across regions and sexes.

  Available options include:

  - `"est_all"`: Estimates separate fixed-effect selectivity parameters
    for each region and sex.

  - `"est_shared_r"`: Shares parameters across regions (sex-specific
    parameters are still estimated).

  - `"est_shared_s"`: Shares parameters across sexes (region-specific
    parameters are still estimated).

  - `"est_shared_r_s"`: Shares parameters across both regions and sexes,
    estimating a single set of fixed-effect parameters.

  - `"est_shared_f_x"`: Shares fixed-effect selectivity parameters with
    another fleet, where `x` is the fleet number to share with. This
    option forces multiple fleets to have identical selectivity curves
    by using the same underlying parameters (e.g., same a50, k, amax
    values). For example, `"est_shared_f_2"` means the current fleet
    will use the same fixed-effect selectivity parameters as fleet 2.
    The reference fleet (fleet x) must use one of the other sharing
    options and cannot itself be sharing with another fleet.

  - `"fix"`: Fixes all selectivity parameters to their initial values
    (no estimation).

  - `"none"`: No selectivity parameters are estimated (equivalent to
    `"fix"`).

- fish_q_spec:

  Character string specifying the structure of fishery catchability
  (`q`) estimation across regions. This controls whether separate or
  shared parameters are used.

  Available options include:

  - `"est_all"`: Estimates separate catchability parameters for each
    region.

  - `"est_shared_r"`: Estimates a single catchability parameter shared
    across all regions.

- fish_sel_devs_spec:

  Character string specifying the structure of process error deviations
  in time-varying fishery selectivity dimensioned by the number of
  fishery fleets. This determines how deviations are estimated across
  regions and sexes.

  Available options include:

  - `"est_all"`: Estimates a separate deviation time series for each
    region and sex.

  - `"est_shared_r"`: Shares deviations across regions (sex-specific
    deviations are still estimated).

  - `"est_shared_s"`: Shares deviations across sexes (region-specific
    deviations are still estimated).

  - `"est_shared_r_s"`: Shares deviations across both regions and sexes,
    estimating a single deviation time series.

  - `"est_shared_a"`: Shares deviations across age blocks.

  - `"est_shared_r_a"`: Shares deviations across regions and age shared
    blocks.

  - `"est_shared_a_s"`: Shares deviations across age shared blocks and
    sexes.

  - `"est_shared_r_a_s"`: Shares deviations across regions, age shared
    blocks, and sexes.

  - `"est_shared_f_x"`: Shares deviations with another fleet, where `x`
    is the fleet number to share with. This option allows multiple
    fleets to use identical deviation parameters, reducing the number of
    parameters to estimate. For example, `"est_shared_f_2"` means the
    current fleet will use the same deviation parameters as fleet 2. The
    reference fleet (fleet x) must use one of the other sharing options
    (`"est_all"`, `"est_shared_r"`, `"est_shared_s"`, or
    `"est_shared_r_s"`) and cannot itself be sharing with another fleet.

  - `"fix"`: Fixes all deviation parameters to zero (no time-variation).

  - `"none"`: No deviation parameters are estimated (equivalent to
    `"fix"`).

  This argument is only used when a continuous time-varying selectivity
  form is specified (e.g., via `cont_tv_fish_sel`).

- corr_opt_semipar:

  Character string specifying which correlation structures to suppress
  when using semi-parametric time-varying selectivity models. Only used
  if `cont_tv_sel` is set to one of `"3dmarg"`, `"3dcond"`, or
  `"2dar1"`.

  This option allows users to turn off estimation of specific
  correlation components in the time-varying selectivity model. This can
  improve stability or enforce assumptions about independence in the
  temporal or age structure.

  Available options:

  - `"corr_zero_y"`: Sets year (temporal) correlations to 0.

  - `"corr_zero_b"`: Sets age correlations to 0.

  - `"corr_zero_y_b"`: Sets both year and bin correlations to 0.

  - `"corr_zero_c"`: Sets cohort correlations to 0. Only valid for
    `cont_tv_sel` = `"3dmarg"` or `"3dcond"`.

  - `"corr_zero_y_c"`: Sets year and cohort correlations to 0. Only
    valid for `cont_tv_sel` = `"3dmarg"` or `"3dcond"`.

  - `"corr_zero_b_c"`: Sets bin (age) and cohort correlations to 0. Only
    valid for `cont_tv_sel` = `"3dmarg"` or `"3dcond"`.

  - `"corr_zero_y_b_c"`: Sets all correlations (year, bin (age), and
    cohort) to 0. Only valid for `cont_tv_sel` = `"3dmarg"` or
    `"3dcond"`; equivalent to an iid structure.

  These correlation-suppression flags are ignored when `cont_tv_sel` is
  set to any other value.

- Use_fish_selex_prior:

  Integer (0 or 1). Flag to enable/disable fishery selectivity priors.
  When set to 1, applies log-normal priors to fishery selectivity
  parameters as specified in `fish_selex_prior`. When set to 0, no
  priors are applied.

- fish_selex_prior:

  Data frame containing prior specifications for fishery selectivity
  parameters. Must include columns: `region` (region index), `fleet`
  (fleet index), `block` (time block index), `sex` (sex index), `par`
  (parameter index), `mu` (prior mean on natural scale), and `sd` (prior
  standard deviation on log scale). Each row specifies a log-normal
  prior N(log(mu), sd) for one selectivity parameter. Only parameters
  with rows in this data frame will have priors applied.

- cont_tv_fish_sel_penalty:

  Whether or not continuous fishery time varying selectivity penalties
  are applied (if cont_tv_fish_sel \> 0)

- fishsel_devs_shared_ages:

  List object for specifying which ages are shared when selectivity
  deviations are semi-parametric (e.g., list(1:5, 6:10, 11:30) specifies
  that ages 1-5, 6-10, and 11-30 have the same deviations.)

- ...:

  Additional arguments specifying starting values for fishery
  selectivity and catchability parameters (fishsel_pe_pars,
  ln_fishsel_devs, ln_fish_fixed_sel_pars, ln_fish_q, fish_q_coeff)

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
