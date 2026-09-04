# Set up movement model inputs and parameter structures

Configures all aspects of spatial movement for the estimation model,
supporting both unstructured Markov transition (`move_type = 0`) and
Continuous Time Markov Chain (`move_type = 1`) formulations, with
optional continuous iid deviations on the movement surface. Validates
all inputs, initializes parameter arrays, constructs TMB/RTMB factor
maps, and populates `input_list$data` accordingly. Must be called after
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

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
  Movement_popblk_spec = "constant",
  Movement_ageblk_spec = "constant",
  Movement_yearblk_spec = "constant",
  Movement_seasblk_spec = "constant",
  Movement_sexblk_spec = "constant",
  cont_vary_movement = "none",
  Movement_cont_pe_pars_spec = "none",
  ctmc_move_dat = NULL,
  adjacency_mat = NULL,
  area_r = rep(1, input_list$data$n_regions),
  diffusion_formula = NULL,
  preference_formula = NULL,
  ctmc_diffusion_bounds = 0,
  ctmc_diffusion_eps = 0.1,
  move_timing = 0,
  ctmc_scale_by_seasdur = 1,
  move_expm_nsub = 0,
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions.

- move_type:

  Integer. Movement model formulation: `0` = unstructured Markov; `1` =
  CTMC. Default `0`.

- do_recruits_move:

  Integer flag. `0` = age-1 fish do not move (default); movement
  deviations and CTMC rows for the minimum age are fixed at zero. `1` =
  recruits participate in movement.

- use_fixed_movement:

  Integer flag. `0` = estimate movement (default); `1` = fix movement
  rates to `Fixed_Movement` and map all movement parameters to `NA`.

- Fixed_Movement:

  Numeric array of externally supplied movement probability matrices,
  dimensioned
  `[n_pop × n_regions × n_regions × n_years × n_seas × n_ages × n_sexes]`.
  Each `[n_regions × n_regions]` slice must be row-stochastic (rows sum
  to 1). Required when `use_fixed_movement = 1`. If `NA` (default), an
  identity matrix (no movement) is constructed internally.

- Use_Movement_Prior:

  Integer flag. `1` = apply Dirichlet priors to movement row
  probabilities; `0` = no priors (default). Requires `Movement_prior`.

- Movement_prior:

  Data frame of Dirichlet prior concentration parameters. Required
  columns: `pop`, `region_from`, `year`, `seas`, `age`, `sex`, and
  `alpha`, where `alpha` is a list-column with each element a numeric
  vector of length `n_regions` giving the Dirichlet concentration for
  transitions out of `region_from`. Values near 1 are uninformative;
  larger values concentrate the prior toward equal movement. Only used
  when `Use_Movement_Prior = 1`.

- Movement_popblk_spec:

  `"constant"` (default, shared across all populations) or a list of
  integer vectors partitioning populations into blocks. Example:
  `list(c(1, 2), 3)` shares parameters for populations 1 and 2 and
  estimates a separate parameter for population 3. Ignored when
  `move_type = 1`.

- Movement_ageblk_spec:

  `"constant"` (default) or a list of integer vectors defining age
  blocks. Example: `list(1:4, 5:10)` creates a juvenile block (ages 1-4)
  and an adult block (ages 5-10). Ignored when `move_type = 1`.

- Movement_yearblk_spec:

  `"constant"` (default) or a list of integer vectors defining year
  blocks for discrete structural breaks in movement. For residual annual
  variation, use `cont_vary_movement` instead. Ignored when
  `move_type = 1`.

- Movement_seasblk_spec:

  `"constant"` (default) or a list of integer vectors defining season
  blocks. Example: `list(c(1, 2), c(3, 4))` groups winter/spring and
  summer/fall. Ignored when `move_type = 1`.

- Movement_sexblk_spec:

  `"constant"` (default, sex-invariant) or a list of integer vectors
  defining sex blocks. Example: `list(1, 2)` estimates sex-specific
  movement independently. Ignored when `move_type = 1`.

- cont_vary_movement:

  Character string specifying the structure of continuous iid movement
  deviations added on top of the fixed-effect movement surface. Default
  `"none"`. Options:

  `"none"`

  :   No deviations.

  `"iid_y"`

  :   Year-varying; shared across pop, age, sex, season.

  `"iid_a"`

  :   Age-varying; shared across pop, year, sex, season.

  `"iid_y_a"`

  :   Year \\\times\\ age.

  `"iid_y_a_s"`

  :   Year \\\times\\ age \\\times\\ sex.

  `"iid_y_seas_a_s"`

  :   Year \\\times\\ season \\\times\\ age \\\times\\ sex.

  `"iid_p_y"`, `"iid_p_a"`, `"iid_p_y_a"`, `"iid_p_y_a_s"`, `"iid_p_y_seas_a_s"`

  :   Population-specific analogs of the above.

- Movement_cont_pe_pars_spec:

  Character string specifying estimation of process-error variance for
  `cont_vary_movement` deviations. One of:

  `"none"`

  :   No process-error parameters; use with
      `cont_vary_movement = "none"`.

  `"fix"`

  :   Parameters initialized but not estimated; fixes deviation variance
      at its starting value.

  `"est_shared"`

  :   Single variance estimated, shared across all dimensions.

  `"est_all"`

  :   All variance parameters estimated independently, dimensioned
      `[n_pop × n_regions × n_seas × n_ages × n_sexes]`.

- ctmc_move_dat:

  Data frame required when `move_type = 1`. Each row corresponds to a
  unique pop-region-year-season-age-sex combination. Required columns:
  `pop`, `regions`, `years`, `seas`, `ages`, `sexes`, plus any covariate
  columns referenced in `diffusion_formula` or `preference_formula`.
  Projection years exceeding `n_years` are automatically capped to the
  final estimation year to prevent spline extrapolation.

- adjacency_mat:

  Square numeric matrix `[n_regions × n_regions]` with 1 indicating an
  allowed transition and 0 indicating no direct connection. The diagonal
  must be 0: residency falls out of the generator, and a non-zero
  diagonal leaves the generator columns summing to something other than
  zero, so the movement matrix loses abundance rather than
  redistributing it. A fully connected matrix is `1 - diag(n_regions)`
  (note that `diag(1, n_regions)` is the identity, not an adjacency
  matrix). Required for `move_type = 1`, where it is validated for
  dimension, 0/1 entries, a zero diagonal, and at least one connection.
  For `move_type = 0` a fully connected matrix is constructed
  automatically.

- area_r:

  Numeric vector of length `n_regions` giving the area of each region,
  used to scale CTMC diffusion rates. Required for `move_type = 1`.
  Default: `rep(1, n_regions)`.

- diffusion_formula:

  R `formula` defining the linear predictor for the CTMC diffusion
  (\\\theta\\) component (e.g., `~ bs(depth, df = 4)`). All
  right-hand-side variables must be present in `ctmc_move_dat`. Required
  for `move_type = 1`.

- preference_formula:

  R `formula` defining the linear predictor for the CTMC
  habitat-preference (taxis, \\\gamma\\) component. All variables must
  be present in `ctmc_move_dat`. Required for `move_type = 1`.

- ctmc_diffusion_bounds:

  How the CTMC generator is kept a valid Metzler matrix (non-negative
  off-diagonal entries) when taxis outweighs diffusion. Writing \\d\\
  for the preference gradient along an edge and \\\theta_j\\ for the
  diffusion rate out of its origin, the accepted values are `"none"` (or
  `0`, the default) for \\\theta_j + d\\ unbounded; `"softplus"` (or
  `1`) for a softplus of \\\theta_j + d\\ of width `ctmc_diffusion_eps`;
  `"clamp"` for its hard positive part; `"upwind"` for the discontinuous
  Galerkin (finite volume) flux \\\theta_j + \max(d, 0)\\; `"barker"`
  for \\\theta_j\\\mathrm{logit}^{-1}(d)\\; and `"logsoftplus"` for
  \\\theta_j \log(1 + e^{d})\\. The last two are positive by
  construction and ignore `ctmc_diffusion_eps` entirely: they have no
  floor, no kink and no gradient-dead region, and they leave the
  diffusion parameter interior when the data want one-way flow. See
  `Get_Movement` for the full comparison.

- ctmc_diffusion_eps:

  Positive numeric width of the softplus applied when
  `ctmc_diffusion_bounds = "softplus"` (default `0.1`). An edge where
  taxis exactly cancels diffusion carries `eps * log(2)`, so this sets a
  floor on exchange as well as smoothing the hinge.

- move_timing:

  Integer flag setting how movement and mortality are sequenced within a
  season. `0` = movement then mortality (default, historical SPoRC
  behavior); `1` = mortality then movement; `2` = continuous, with
  movement and mortality acting simultaneously via the matrix
  exponential of \\Q\Delta - \mathrm{diag}(Z)\\. `move_timing = 2`
  requires an estimated CTMC generator, i.e. `move_type = 1` and
  `use_fixed_movement = 0`.

- ctmc_scale_by_seasdur:

  Integer flag controlling the time units of the CTMC generator. `1`
  (default) treats \\Q\\ as an annual rate, exponentiating \\Q \cdot
  \mathrm{seasdur}\[s\]\\ in each season so that movement and mortality
  share time units. `0` exponentiates \\Q\\ once per season regardless
  of duration. Only has an effect when `move_type = 1` and `n_seas > 1`;
  forced to `1` when `move_timing = 2`.

- move_expm_nsub:

  Integer controlling how matrix exponentials of the CTMC generator are
  evaluated, both when converting \\Q\\ to movement fractions and inside
  the `move_timing = 2` seasonal operators. `0` (default) uses
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html). A
  power of two \\n \ge 1\\ instead uses \\n\\ implicit (backward Euler)
  substeps, \\(I - A/n)^{-n}\\, evaluated as one linear solve plus
  \\\log_2 n\\ squarings, which is why \\n\\ must be a power of two. The
  implicit form has a much cheaper reverse-mode adjoint than a matrix
  exponential, so the gradient is several times faster, but it is a
  first-order approximation: \\n = 1\\ is plain `solve(I - A)` and is an
  approximation.

- ...:

  Optional starting value overrides, passed by name. Recognized
  arguments:

  `move_pars`

  :   Array
      `[n_pop × n_regions × (n_regions-1) × n_years × n_seas × n_ages × n_sexes]`.
      Default: `0` (equal movement on logit scale).

  `log_move_diffusion_pars`

  :   Vector of length `n_theta`. Default: `log(0.1)`.

  `move_preference_pars`

  :   Vector of length `n_gamma`. Default: `0`.

  `move_devs`

  :   Array
      `[n_pop × n_regions × (n_regions-1) × (n_years + n_proj_yrs_devs) × n_seas × n_ages × n_sexes]`.
      Default: `0`.

  `move_pe_pars`

  :   Array `[n_pop × n_regions × n_seas × n_ages × n_sexes]`. Default:
      `0`.

## Value

The input `input_list` with `$data`, `$par`, and `$map` updated. Key
additions to `$data` include `move_type`, `use_fixed_movement`,
`Fixed_Movement`, `adjacency_mat`, `adjacency_collapsed`, `area_r`,
`ctmc_move_dat`, `diffusion_formula`, `preference_formula`, and
`cont_vary_movement` (stored as an integer code). Parameter arrays
`move_pars`, `log_move_diffusion_pars`, `move_preference_pars`,
`move_devs`, and `move_pe_pars` are added to `$par`, with corresponding
factor maps in `$map`.

## Unstructured Markov movement (`move_type = 0`)

Transition probabilities from region \\r\\ to all other regions are
parameterized via a multinomial logit with a reference-cell constraint.
The parameter array `move_pars` has dimensions
`[n_pop × n_regions × (n_regions - 1) × n_years × n_seas × n_ages × n_sexes]`.
Block specifications (`Movement_*blk_spec`) control sharing: indices
within the same block receive the same TMB factor level and are
estimated as a single free parameter. A fully connected adjacency matrix
is constructed automatically. Blocked and continuous time-varying
movement can be combined: use `Movement_yearblk_spec` for discrete
structural breaks and `cont_vary_movement` for residual year-to-year
variation.

## CTMC movement (`move_type = 1`)

The instantaneous rate matrix \\Q\\ is decomposed into diffusion
(\\\theta\\) and preference (\\\gamma\\) components following Thorson et
al. Design matrices for both components are derived from
`diffusion_formula` and `preference_formula` evaluated on
`ctmc_move_dat`. The discrete-time movement matrix for each time step is
\\\exp(Q \Delta t)\\. Parameter blocking is not supported for CTMC; all
`Movement_*blk_spec` arguments must remain `"constant"`. Structural
variation across populations, ages, sexes, or seasons should instead be
introduced through formula covariates in `ctmc_move_dat`.

## Continuous movement deviations

IID deviations (`move_devs`) are added to the movement logit surface
(unstructured Markov) or log-rate surface (CTMC) before computing
probabilities. Deviations are penalized as normal random effects; the
variance is optionally estimated via `Movement_cont_pe_pars_spec`. If
`do_recruits_move = 0`, age-1 deviations are fixed at zero.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
