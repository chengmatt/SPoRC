# Setup Movement Processes for SPoRC

Configure movement model components (unstructured Markov or CTMC) and
populate the model input and parameter lists with the appropriate data
structures and starting values.

## Usage

``` r
Setup_Mod_Movement(
  input_list,
  move_type = 0,
  do_recruits_move = 0,
  use_fixed_movement = 0,
  Fixed_Movement = NA,
  Use_Movement_Prior = 0,
  Movement_prior = NULL,
  Movement_ageblk_spec = "constant",
  Movement_yearblk_spec = "constant",
  Movement_sexblk_spec = "constant",
  cont_vary_movement = "none",
  Movement_cont_pe_pars_spec = "none",
  ctmc_move_dat = NULL,
  adjacency_mat = NULL,
  area_r = rep(1, input_list$data$n_regions),
  diffusion_formula = NULL,
  preference_formula = NULL,
  ctmc_diffusion_bounds = 0,
  ...
)
```

## Arguments

- input_list:

  List containing data, parameter, and map lists for the model. This
  object is updated in-place and returned by the function.

- move_type:

  Integer indicating the movement model type:

  - `0` – unstructured Markov (parameters transformed via a multinomial
    logit),

  - `1` – Continuous Time Markov Chain (CTMC) formulation.

  Default is `0`.

- do_recruits_move:

  Integer flag (0 or 1) indicating whether recruits move. Default is 0
  (recruits do not move).

- use_fixed_movement:

  Integer flag (0 or 1) indicating whether to use a fixed movement
  matrix (1) or estimate movement parameters (0). Default is 0.

- Fixed_Movement:

  Numeric array specifying a fixed movement rate/matrix. It must be
  dimensioned as `[n_regions, n_regions, n_years, n_ages, n_sexes]`. If
  `NA` (the default), a neutral array of ones will be created
  internally.

- Use_Movement_Prior:

  Integer flag (0 or 1) indicating whether to use movement priors.
  Default is 0 (priors not used).

- Movement_prior:

  Optional data.frame providing informative priors for movement.
  Required columns are `region_from`, `age`, `sex`, and `alpha`, where
  `alpha` is a list-column and each element is a numeric vector of
  length `n_regions` containing prior concentration parameters for
  movement from the specified region. If `NULL` (default), no movement
  prior is used.

- Movement_ageblk_spec:

  Only applicable for move_type = 0. Either:

  - Character string `"constant"` for age-invariant movement (default),
    or

  - A `list` of integer vectors specifying age blocks that share
    parameters.

  Example: `list(c(1:6), c(7:10), c(11:n_ages))` makes ages 1–6, 7–10,
  and 11–n_ages share parameters. To indicate full age invariance, use
  `"constant"` or `list(c(1:n_ages))`.

- Movement_yearblk_spec:

  Only applicable for move_type = 0. Either:

  - Character string `"constant"` for time-invariant movement (default),
    or

  - A `list` of integer vectors specifying year blocks that share
    movement parameters.

- Movement_sexblk_spec:

  Only applicable for move_type = 0. Either:

  - Character string `"constant"` for sex-invariant movement (default),
    or

  - A `list` of integer vectors specifying sex blocks that share
    movement parameters.

- cont_vary_movement:

  Character string specifying continuous varying movement type.
  Available options:

  - `"none"` (no continuous variation)

  - `"iid_y"` (iid deviations by year)

  - `"iid_a"` (iid deviations by age)

  - `"iid_y_a"` (iid deviations by year and age)

  - `"iid_y_s"` (iid deviations by year and sex)

  - `"iid_a_s"` (iid deviations by age and sex)

  - `"iid_y_a_s"` (iid deviations by year, age, and sex)

  Default is `"none"`.

- Movement_cont_pe_pars_spec:

  Character string specifying how process-error parameters for
  continuous-varying movement are shared or estimated. Available
  options:

  - `"est_shared_r"` – estimate shared across regions,

  - `"est_shared_a"` – estimate shared across ages,

  - `"est_shared_s"` – estimate shared across sexes,

  - `"est_shared_r_a"`, `"est_shared_a_s"`, `"est_shared_r_s"`,

  - `"est_shared_r_a_s"` – combinations of shared structure,

  - `"est_all"` – estimate all process-error parameters independently,

  - `"fix"` – treat process-error parameters as fixed (not estimated),

  - `"none"` – no process-error parameters (default).

  Default is `"none"`.

- ctmc_move_dat:

  Data.frame with CTMC covariates used to build design matrices for
  diffusion and preference. Required columns (when `move_type == 1`)
  include `regions`, `years`, `ages`, and `sexes`, plus any covariates
  referenced in `diffusion_formula` and `preference_formula`. Can
  include projection years (years \> n_yrs) with projected covariate
  values. Year effects in formulas (e.g., splines) are automatically
  capped at `n_yrs`

- adjacency_mat:

  Square adjacency matrix (`n_regions x n_regions`) for CTMC movement
  (used when `move_type == 1`). Non-zero entries indicate allowed
  transitions between region pairs; zero entries indicate transitions
  that are not allowed.

- area_r:

  Numeric vector of region areas of length `n_regions`. Required for
  CTMC movement to convert rates to per-area or per-distance values.
  Defaults to `rep(1, n_regions)`.

- diffusion_formula:

  An R formula describing the linear predictor for diffusion rates in
  the CTMC model (required when `move_type == 1`). Variables used in the
  formula must be present as columns in `ctmc_move_dat`.

- preference_formula:

  An R formula describing the linear predictor for preference in the
  CTMC model (required when `move_type == 1`). Variables used in the
  formula must be present as columns in `ctmc_move_dat`.

- ctmc_diffusion_bounds:

  Numeric indicating whether diffusion bounds are placed to ensure that
  the generator matrix is Metzler. 0 == no bounds (default), 1 == bounds
  enforced, representing slipstream diffusion (i.e., residual preference
  gradients get added to diffusion to ensure Metzler matrix).

- ...:

  Additional named starting values that may be supplied. Typical names:
  `move_pars`, `move_devs`, `move_pe_pars`, `log_move_diffusion_pars`,
  `move_preference_pars`, etc. If not supplied, sensible defaults are
  created.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
