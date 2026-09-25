# Set up movement model inputs and parameter structures

Sets up unstructured Markov transition movement (`move_type = 0`) or a
continuous time Markov chain (`move_type = 1`), with optional iid
deviations on the movement surface, and builds the parameter arrays and
factor maps. Call after
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

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- move_type:

  Integer. `0` (default) unstructured Markov, `1` CTMC.

- do_recruits_move:

  Integer flag. `0` (default) fixes the movement deviations and CTMC
  rows at the minimum age to zero, `1` moves recruits.

- use_fixed_movement:

  Integer flag. `0` (default) estimates movement, `1` fixes it at
  `Fixed_Movement` and maps every movement parameter to `NA`.

- Fixed_Movement:

  Movement probability array
  `[n_pop × n_regions × n_regions × n_years × n_seas × n_ages × n_sexes]`,
  each `[n_regions × n_regions]` slice row-stochastic. Required when
  `use_fixed_movement = 1`. `NA` (default) builds an identity matrix.

- Use_Movement_Prior:

  Integer flag, `1` for Dirichlet priors on the movement rows. Default
  `0`.

- Movement_prior:

  Data frame with columns `pop`, `region_from`, `year`, `seas`, `age`,
  `sex` and `alpha`, the last a list-column of length-`n_regions`
  concentrations for transitions out of `region_from`. Values near 1 are
  uninformative, larger ones concentrate toward equal movement. Read
  when `Use_Movement_Prior = 1`.

- Movement_popblk_spec, Movement_ageblk_spec, Movement_yearblk_spec,
  Movement_seasblk_spec, Movement_sexblk_spec:

  Blocking across populations, ages, years, seasons and sexes:
  `"constant"` (default) or a list of integer vectors, e.g.
  `list(c(1, 2), 3)` for populations, `list(1:4, 5:10)` for a juvenile
  and an adult block, or `list(1, 2)` for sex-specific movement. Use
  `Movement_yearblk_spec` for structural breaks and `cont_vary_movement`
  for residual annual variation. All are ignored when `move_type = 1`.

- cont_vary_movement:

  Structure of the continuous deviations on the fixed-effect movement
  surface. `"none"` (default), or `"iid_"` followed by the dims they
  vary over, any of p (population), y (year), seas (season), a (age) and
  s (sex) in any order: `"iid_y"` is one deviation per year and region
  pair, or per year and region under CTMC movement, shared across
  everything else, and `"iid_p_y_seas_a_s"` varies by every dim. A dim
  left out shares one deviation across it. They are random effects with
  `Movement_cont_pe_pars_spec` estimating the sd and
  `random = "move_devs"` in
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).
  `"dsem"` instead hands their density to the arrows given to
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md),
  one series per origin and destination (per region under CTMC, whose
  deviations hold no destination) and per level of every other dim with
  more than one, which it names itself since a deviation shared across a
  dim cannot be linked; `move_pe_pars` are then read by nothing.

- Movement_cont_pe_pars_spec:

  Estimation of the process error variance for the `cont_vary_movement`
  deviations. `"none"` creates no parameters and pairs with
  `cont_vary_movement = "none"`, `"fix"` holds the variance at its
  starting value, `"est_shared"` estimates one shared value, and
  `"est_all"` estimates
  `[n_pop × n_regions × n_seas × n_ages × n_sexes]` independently.

- ctmc_move_dat:

  Data frame required when `move_type = 1`, one row per population,
  region, year, season, age and sex, with columns `pop`, `regions`,
  `years`, `seas`, `ages`, `sexes` and any covariates the formulas name.
  Projection years beyond `n_years` are capped at the final estimation
  year to prevent spline extrapolation.

- adjacency_mat:

  Square `[n_regions × n_regions]` matrix, 1 for an allowed transition
  and 0 for none. The diagonal must be 0: residency falls out of the
  generator, and a non-zero diagonal leaves the generator columns
  summing to something other than zero, so the movement matrix loses
  abundance rather than redistributing it. A fully connected matrix is
  `1 - diag(n_regions)` (`diag(1, n_regions)` is the identity, not an
  adjacency matrix). Required under `move_type = 1`, where it is
  validated for dimension, 0/1 entries, a zero diagonal and at least one
  connection; built automatically under `move_type = 0`.

- area_r:

  Numeric vector `[n_regions]` of region areas, used to scale the CTMC
  diffusion rates. Required under `move_type = 1`. Default
  `rep(1, n_regions)`.

- diffusion_formula:

  Formula for the CTMC diffusion (\\\theta\\) linear predictor, e.g.
  `~ bs(depth, df = 4)`. Every right-hand-side variable must be in
  `ctmc_move_dat`. Required under `move_type = 1`.

- preference_formula:

  Formula for the CTMC preference (taxis, \\\gamma\\) linear predictor,
  on the same terms. Required under `move_type = 1`.

- ctmc_diffusion_bounds:

  How the CTMC generator is kept a valid Metzler matrix when taxis
  outweighs diffusion. `"softplus"` takes a softplus of \\\theta_j + d\\
  of width `ctmc_diffusion_eps`; `"upwind"` takes the finite volume flux
  \\\theta_j + \max(d, 0)\\, which keeps diffusion whole and adds only
  the down-gradient taxis, so positivity never depends on the two
  cancelling.

- ctmc_diffusion_eps:

  Positive width of the softplus under
  `ctmc_diffusion_bounds = "softplus"`. Default `0.1`. An edge where
  taxis exactly cancels diffusion has `eps * log(2)`, so this is a floor
  on exchange as well as a smoothing constant.

- move_timing:

  How movement and mortality are sequenced within a season. `0`
  (default) moves then kills, `1` kills then moves, and `2` runs the two
  together through the matrix exponential of \\Q\Delta -
  \mathrm{diag}(Z)\\. `2` needs an estimated CTMC generator, so
  `move_type = 1` and `use_fixed_movement = 0`.

- ctmc_scale_by_seasdur:

  Integer flag for the time units of the CTMC generator. `1` (default)
  treats \\Q\\ as an annual rate and exponentiates \\Q \cdot
  \mathrm{seasdur}\[s\]\\ each season, so movement and mortality share
  time units; `0` exponentiates \\Q\\ once per season whatever its
  duration. Only matters under `move_type = 1` with `n_seas > 1`, and is
  forced to `1` under `move_timing = 2`.

- move_expm_nsub:

  How matrix exponentials of the generator are evaluated, both
  converting \\Q\\ to movement fractions and inside the
  `move_timing = 2` seasonal operators. `0` (default) uses
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html). A
  power of two \\n \ge 1\\ uses \\n\\ implicit backward Euler substeps,
  \\(I - A/n)^{-n}\\, as one linear solve plus \\\log_2 n\\ squarings,
  which is why \\n\\ must be a power of two. Its reverse-mode derivative
  is much cheaper, so the gradient is several times faster, but it is a
  first-order approximation and \\n = 1\\ is plain `solve(I - A)`.

- ...:

  Optional starting values by name: `move_pars`
  `[n_pop × n_regions × (n_regions-1) × n_years × n_seas × n_ages × n_sexes]`,
  default `0`; `log_move_diffusion_pars` of length `n_theta`, default
  `log(0.1)`; `move_preference_pars` of length `n_gamma`, default `0`;
  `move_devs`, shaped as `move_pars` with `n_years + n_proj_yrs_devs`
  years and the third dim `1` under `move_type = 1`, default `0`; and
  `move_pe_pars` `[n_pop × n_regions × n_seas × n_ages × n_sexes]`,
  default `0`.

## Value

`input_list` with `$data`, `$par` and `$map` updated. `$data` gains
`move_type`, `use_fixed_movement`, `Fixed_Movement`, `adjacency_mat`,
`adjacency_collapsed`, `area_r`, `ctmc_move_dat`, `diffusion_formula`,
`preference_formula` and `cont_vary_movement` as its form string.
`move_pars`, `log_move_diffusion_pars`, `move_preference_pars`,
`move_devs` and `move_pe_pars` go into `$par` with their factor maps in
`$map`.

## Unstructured Markov movement (`move_type = 0`)

Transitions out of region \\r\\ are a multinomial logit with a reference
cell, so `move_pars` is
`[n_pop × n_regions × (n_regions - 1) × n_years × n_seas × n_ages × n_sexes]`.
The `Movement_*blk_spec` arguments share parameters: indices in one
block take the same factor level. A fully connected adjacency matrix is
built automatically. Blocks and continuous time variation combine: use
`Movement_yearblk_spec` for structural breaks and `cont_vary_movement`
for residual year-to-year variation.

## CTMC movement (`move_type = 1`)

The rate matrix \\Q\\ is decomposed into diffusion (\\\theta\\) and
preference (\\\gamma\\), with design matrices from `diffusion_formula`
and `preference_formula` evaluated on `ctmc_move_dat`, and each time
step's movement matrix is \\\exp(Q \Delta t)\\. Blocking is not
supported, so every `Movement_*blk_spec` must stay `"constant"`; put
structure across populations, ages, sexes or seasons into formula
covariates instead.

## Continuous movement deviations

Deviations are added to the movement logit surface under unstructured
movement, or to each region's preference under CTMC, before
probabilities are computed, and are penalized as normal random effects
whose variance `Movement_cont_pe_pars_spec` can estimate. Age-1
deviations are fixed at zero when `do_recruits_move = 0`.

The two types size the deviations differently. Unstructured movement
holds one per origin and destination pair,
`[n_regions x (n_regions - 1)]`, while the CTMC holds one per region,
`[n_regions x 1]`, since preference is a surface over regions rather
than a rate along an edge. A CTMC deviation raises the rates into its
region and lowers those out of it on every edge at once. Only
differences in preference reach the generator, so a constant added to
every region's deviation leaves movement unchanged and nothing but the
process error penalty holds the level of the field down.

Because the deviations are preference they enter the generator
additively, and one larger than an edge's diffusion rate drives that
rate negative. Under `ctmc_diffusion_bounds = "none"` the movement
fractions then leave `[0, 1]` while still summing to one, so use
`"upwind"` or `"softplus"` whenever the deviations are estimated.

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
