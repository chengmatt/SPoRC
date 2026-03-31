# Construct Movement Matrices for Unstructured or CTMC Movement

Generates movement matrices for a population model based on either
unstructured multinomial logit movement or a Continuous Time Markov
Chain (CTMC) formulation. Also calculates a movement penalty if
applicable. For CTMC movement, projection years are supported: base
parameters (preference/diffusion) can either be frozen at the last
historical year or extended via user-provided covariates in
ctmc_move_dat.

## Usage

``` r
Get_Movement(
  move_type,
  do_recruits_move,
  n_regions,
  n_yrs,
  n_proj_yrs_devs,
  n_ages,
  n_sexes,
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
  ctmc_diffusion_bounds
)
```

## Arguments

- move_type:

  Integer flag indicating movement type: 0 = unstructured Markov, 1 =
  CTMC movement.

- do_recruits_move:

  Integer flag: 0 = recruits do not move, 1 = recruits move.

- n_regions:

  Number of spatial regions.

- n_yrs:

  Number of years in the observed data.

- n_proj_yrs_devs:

  Number of projected years for deviations.

- n_ages:

  Number of age classes.

- n_sexes:

  Number of sexes.

- move_pars:

  Array of movement parameters for unstructured movement.

- move_devs:

  Array of movement deviations (applies to both unstructured and CTMC
  movement).

- use_fixed_movement:

  Integer flag: 0 = estimate movement, 1 = use fixed matrix.

- Fixed_Movement:

  Optional fixed movement matrix.

- ctmc_move_dat:

  Data.frame with CTMC covariates used to build design matrices for
  diffusion and preference. Required columns (when `move_type == 1`)
  include `regions`, `years`, `ages`, and `sexes`, plus any covariates
  referenced in `diffusion_formula` and `preference_formula`. Can
  include projection years (years \> n_yrs) with projected covariate
  values. Year effects in formulas (e.g., splines) are automatically
  capped at `n_yrs` to prevent extrapolation, while other covariates use
  their actual projected values.

- preference_formula:

  R formula specifying preference covariates for CTMC movement.

- diffusion_formula:

  R formula specifying diffusion covariates for CTMC movement.

- log_move_diffusion_pars:

  Log-transformed diffusion parameters for CTMC movement.

- move_preference_pars:

  Preference parameters for CTMC movement.

- area_r:

  Vector of areas for each region (used for scaling diffusion rates).

- adjacency_mat:

  Square adjacency matrix defining connectivity between regions for CTMC
  movement.

- ctmc_diffusion_bounds:

  Integer flag: 1 = apply diffusion bounds to generator matrix, 0 = no
  bounds.

## Value

A list with components:

- `Movement`:

  Array of movement fractions for each stratum (from regions × to
  regions × years × ages × sexes).

- `Mrate`:

  Instantaneous movement rate matrix if CTMC movement is used (from
  regions × to regions × years × ages × sexes); otherwise NULL.

- `move_pen`:

  Numeric value of movement penalty calculated from preference
  parameters (for CTMC only).
