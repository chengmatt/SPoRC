# Construct Movement Matrices for Unstructured or CTMC Movement

Generates movement matrices for a population model based on either
unstructured multinomial logit movement or a Continuous Time Markov
Chain (CTMC) formulation. Also calculates a movement penalty if
applicable. For CTMC movement, projection years are supported: covariate
lookups are capped at the last historical year (`n_yrs`), so CTMC
parameters are frozen at their final historical values during projection
unless `ctmc_move_dat` is extended with projection-year rows.

## Usage

``` r
Get_Movement(
  move_type,
  do_recruits_move,
  n_pop,
  n_regions,
  n_yrs,
  n_proj_yrs_devs,
  n_ages,
  n_sexes,
  n_seas,
  move_pars,
  move_devs,
  use_fixed_movement,
  Fixed_Movement = NULL,
  ctmc_move_dat = NULL,
  preference_formula = NULL,
  diffusion_formula = NULL,
  log_move_diffusion_pars,
  move_preference_pars,
  area_r,
  adjacency_mat,
  ctmc_diffusion_bounds,
  seasdur = rep(1, n_seas),
  ctmc_scale_by_seasdur = 0
)
```

## Arguments

- move_type:

  Integer flag indicating movement type: 0 = unstructured Markov, 1 =
  CTMC movement.

- do_recruits_move:

  Integer flag: 0 = recruits (age 1) do not move, 1 = recruits move.

- n_pop:

  Number of populations.

- n_regions:

  Number of spatial regions.

- n_yrs:

  Number of years in the observed (historical) data.

- n_proj_yrs_devs:

  Number of projected years. Extends the year dimension of the movement
  array beyond `n_yrs`; CTMC covariate lookups are capped at `n_yrs`
  unless `ctmc_move_dat` contains projection-year rows.

- n_ages:

  Number of age classes.

- n_sexes:

  Number of sexes.

- n_seas:

  Number of seasons.

- move_pars:

  Array of movement parameters for unstructured (multinomial logit)
  movement, dimensioned
  `[pop, from_region, counter, year, seas, age, sex]`. `counter` indexes
  the `n_regions - 1` non-reference destinations. Ignored when
  `move_type == 1` or `use_fixed_movement == 1`.

- move_devs:

  Array of origin-destination movement deviations, dimensioned
  `[pop, origin_region, counter, year, seas, age, sex]`, where `counter`
  indexes adjacent non-diagonal destinations from each origin region.
  Applied as multiplicative offsets (`exp(move_devs)`) to off-diagonal
  diffusion rates in CTMC movement, or as additive offsets to
  logit-scale parameters in unstructured movement. Always indexed on the
  actual (possibly projected) year.

- use_fixed_movement:

  Integer flag: 0 = estimate movement, 1 = use fixed matrix.

- Fixed_Movement:

  Fixed movement matrix used when `use_fixed_movement == 1`, dimensioned
  `[pop, from_region, to_region, year, seas, age, sex]`. Ignored
  otherwise.

- ctmc_move_dat:

  Data.frame of CTMC covariates. Required when `move_type == 1`. Must
  include columns `pop`, `regions`, `years`, `seas`, `ages`, `sexes`,
  and any covariates referenced in `diffusion_formula` or
  `preference_formula`.

- preference_formula:

  R formula specifying preference (taxis) covariates for CTMC movement.
  Required when `move_type == 1`.

- diffusion_formula:

  R formula specifying diffusion covariates for CTMC movement. Required
  when `move_type == 1`.

- log_move_diffusion_pars:

  Log-scale diffusion parameters (\\\theta_k\\) for CTMC movement.
  Exponentiated and squared internally (`exp(2 * log_theta)`). Required
  when `move_type == 1`.

- move_preference_pars:

  Preference (taxis) parameters (\\\gamma_k\\) for CTMC movement, on the
  natural scale. Required when `move_type == 1`.

- area_r:

  Numeric vector of region areas (length `n_regions`) used to scale
  diffusion rates. Required when `move_type == 1`.

- adjacency_mat:

  Square `[n_regions x n_regions]` matrix defining connectivity among
  regions (1 = adjacent, 0 = not adjacent). Required when
  `move_type == 1`.

- ctmc_diffusion_bounds:

  Integer flag: 1 = shift diffusion columns to ensure all off-diagonal
  generator matrix entries are non-negative (valid generator); 0 = no
  bounds applied.

- seasdur:

  Numeric vector of length `n_seas` giving season durations (summing to
  1). Used to scale the CTMC generator when
  `ctmc_scale_by_seasdur == 1`. Defaults to `rep(1, n_seas)`, which
  reproduces the unscaled behaviour.

- ctmc_scale_by_seasdur:

  Integer flag controlling the time units of the CTMC generator. `1` =
  treat \\Q\\ as an annual rate and exponentiate \\Q \cdot
  \mathrm{seasdur}\[s\]\\ for each season; `0` = treat \\Q\\ as a per
  season rate and exponentiate it once per season irrespective of season
  duration. Only affects `move_type == 1` with `n_seas > 1`. Defaults to
  `0` here so that callers passing an unscaled generator get the
  arithmetic they expect; the user facing default is `1`, set by
  `Setup_Mod_Movement`.

  Whether this flag changes the fit or only reparameterises it depends
  on whether the generator varies by season. With a season-agnostic
  generator it matters: the seasonal steps commute, so scaling on
  composes across the year to \\\exp(Q)\\ regardless of `n_seas`,
  whereas scaling off composes to \\\exp(n\_{seas} Q)\\ and inflates
  movement as the seasonal time step shrinks. With a season-varying
  generator the scaling is absorbed — \\Q\\ is linear in \\\theta\\, so
  scaling \\Q\\ by `seasdur[s]` equals shifting
  `log_move_diffusion_pars` by \\\tfrac{1}{2}\log
  \mathrm{seasdur}\[s\]\\ — but only if *both* formulas carry a season
  term, since the flag scales diffusion and taxis together while only
  \\\theta\\ can absorb it.

  Note that even where the scaling is absorbed, it changes what the
  estimated diffusion means (per-season step vs annual rate, hence not
  comparable across seasons of unequal length), and the Dirichlet
  movement prior is evaluated on annual fractions \\\exp(Q)\\
  irrespective of this flag. Under `move_timing = 2` the flag governs
  only the reported `Movement` diagnostic: the dynamics take `Mrate` and
  apply `seasdur[seas]` themselves.

## Value

A list with components:

- `Movement`:

  Array of movement fractions dimensioned
  `[pop, from_region, to_region, year, seas, age, sex]`. Populated for
  all three movement configurations.

- `Mrate`:

  Instantaneous rate matrix (generator \\Q\\) dimensioned
  `[pop, from_region, to_region, year, seas, age, sex]`. Populated only
  when `move_type == 1` and `use_fixed_movement == 0`; `NULL` otherwise.
  Stored *unscaled* by season duration, so that consumers combining it
  with mortality must apply `seasdur[seas]` themselves (see
  `build_seas_operator`).

- `move_pen`:

  Numeric movement penalty used for regularization. For CTMC movement,
  equal to \\\sum\_{\text{strata}} (\sum_r \gamma\_{z,r})^2\\, which
  penalizes large net preference values across regions. Zero for
  unstructured or fixed movement.

## Details

Three movement configurations are supported:

1.  Fixed movement (`use_fixed_movement == 1`): `Fixed_Movement` is used
    directly with no estimation.

2.  Unstructured multinomial logit movement (`move_type == 0`): movement
    fractions are estimated via a softmax transform of
    `move_pars + move_devs`.

3.  CTMC movement (`move_type == 1`): a generator matrix \\Q = D + Z\\
    is constructed from diffusion (\\D\\) and taxis (\\Z\\) components,
    then exponentiated via `Matrix::expm(Q)` to obtain movement
    fractions. During projection years (`y > n_yrs`), covariate lookups
    are capped at `n_yrs` (i.e., CTMC base parameters are frozen at
    their last historical values), while `move_devs` continue to use the
    actual projected year index.
