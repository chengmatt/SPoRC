# Construct Movement Matrices for Unstructured or CTMC Movement

Movement fractions under unstructured multinomial logit movement, a CTMC
generator, or a fixed matrix, plus the movement penalty. Under CTMC
movement the covariate lookups are capped at `n_yrs`, so the parameters
are frozen at their last historical values through the projection unless
`ctmc_move_dat` holds projection-year rows.

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
  ctmc_diffusion_eps = 0.1,
  seasdur = rep(1, n_seas),
  ctmc_scale_by_seasdur = 0,
  expm_nsub = 0
)
```

## Arguments

- move_type:

  Integer: 0 unstructured Markov, 1 CTMC.

- do_recruits_move:

  Integer: 0 keeps recruits (age 1) in place, 1 moves them.

- n_regions, n_ages, n_sexes, n_seas, n_pop:

  Model dimensions.

- n_yrs:

  Number of observed (historical) years.

- n_proj_yrs_devs:

  Number of projected years, extending the year dim of the movement
  array beyond `n_yrs`. CTMC covariate lookups stay capped at `n_yrs`
  unless `ctmc_move_dat` holds projection-year rows.

- move_pars:

  Unstructured movement parameters
  `[pop, from_region, counter, year, seas, age, sex]`, `counter`
  indexing the `n_regions - 1` non-reference destinations. Ignored when
  `move_type == 1` or `use_fixed_movement == 1`.

- move_devs:

  Movement deviations `[pop, region, counter, year, seas, age, sex]`,
  always indexed on the actual (possibly projected) year. Under
  unstructured movement `counter` indexes the non-reference destinations
  and the deviations are additive on the logit scale. Under CTMC
  movement they sit on each region's preference instead, so `counter`
  has length one and `region` is the region whose preference is shifted.

- use_fixed_movement:

  Integer: 0 estimates movement, 1 uses `Fixed_Movement`.

- Fixed_Movement:

  Fixed movement matrix
  `[pop, from_region, to_region, year, seas, age, sex]`, read when
  `use_fixed_movement == 1`.

- ctmc_move_dat:

  Data frame of CTMC covariates, required when `move_type == 1`, with
  columns `pop`, `regions`, `years`, `seas`, `ages`, `sexes` and any
  covariate the formulas name.

- preference_formula, diffusion_formula:

  Formulas for the preference (taxis) and diffusion covariates. Required
  when `move_type == 1`.

- log_move_diffusion_pars:

  Log-scale diffusion parameters (\\\theta_k\\), exponentiated and
  squared internally as `exp(2 * log_theta)`. Required when
  `move_type == 1`.

- move_preference_pars:

  Preference (taxis) parameters (\\\gamma_k\\) on the natural scale.
  Required when `move_type == 1`.

- area_r:

  Numeric vector `[n_regions]` of region areas, scaling the diffusion
  rates. Required when `move_type == 1`.

- adjacency_mat:

  Square `[n_regions x n_regions]` connectivity matrix, 1 for adjacent
  and 0 for not. Required when `move_type == 1`.

- ctmc_diffusion_bounds:

  How the off-diagonal generator entries are kept non-negative. Every
  form is evaluated on the adjacency edges only, so non-edges stay
  exactly zero. With \\d\\ the preference gradient \\\gamma_i -
  \gamma_j\\ along the edge from \\j\\ to \\i\\ and \\\theta_j\\ the
  diffusion rate out of \\j\\: `"none"` (or `0`) is \\q = \theta_j +
  d\\, valid only where diffusion outweighs taxis everywhere;
  `"softplus"` (or `1`) is a softplus of that of width
  `ctmc_diffusion_eps`, smooth, but an edge where taxis cancels
  diffusion has a floor of `eps * log(2)`, so the width is a minimum
  exchange rate and not only a smoothing constant; `"upwind"` (or `2`)
  is the finite volume upwind flux \\q = \theta_j + \max(d, 0)\\, which
  keeps diffusion whole and adds only the down-gradient half of the
  taxis flux, so positivity never depends on the two cancelling.

- ctmc_diffusion_eps:

  Positive width of the softplus under
  `ctmc_diffusion_bounds = "softplus"`. Default 0.1.

- seasdur:

  Numeric vector `[n_seas]` of season durations summing to 1, used to
  scale the generator when `ctmc_scale_by_seasdur == 1`. Defaults to
  `rep(1, n_seas)`, the unscaled behavior.

- ctmc_scale_by_seasdur:

  Integer flag for the generator's time units. `1` treats \\Q\\ as an
  annual rate and exponentiates \\Q \cdot \mathrm{seasdur}\[s\]\\ each
  season; `0` treats it as a per season rate and exponentiates it once
  per season. Only matters under `move_type == 1` with `n_seas > 1`.
  Defaults to `0` here so a caller passing an unscaled generator gets
  the arithmetic it expects; the user facing default is `1`, set by
  `Setup_Mod_Movement`.

- expm_nsub:

  How the generator is exponentiated into movement fractions. `0`
  (default) uses
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html);
  \\n \ge 1\\ uses the implicit backward Euler scheme \\(I -
  Q\Delta/n)^{-n}\\, cheaper to differentiate but first order. Read when
  `move_type == 1` and `use_fixed_movement == 0`. See
  [`mat_exp`](https://chengmatt.github.io/SPoRC/dev/reference/mat_exp.md).

## Value

A list with `Movement`, the movement fractions
`[pop, from_region, to_region, year, seas, age, sex]`, populated under
every configuration; `Mrate`, the generator \\Q\\ on the same dims,
populated only when `move_type == 1` and `use_fixed_movement == 0` and
stored unscaled by season duration, so a consumer combining it with
mortality applies `seasdur[seas]` itself (see `build_seas_operator`);
and `move_pen`, the penalty, which under CTMC movement is \\\sum_k
\gamma_k^2\\, a ridge on the preference coefficients applied once rather
than per stratum that pins the otherwise unidentified level and spread,
and is zero for unstructured or fixed movement.

## Details

Fixed movement (`use_fixed_movement == 1`) uses `Fixed_Movement`
directly. Unstructured movement (`move_type == 0`) takes a softmax of
`move_pars + move_devs`. CTMC movement (`move_type == 1`) builds a
generator \\Q = D + Z\\ from its diffusion and taxis components and
exponentiates it through
[`mat_exp`](https://chengmatt.github.io/SPoRC/dev/reference/mat_exp.md).

Under CTMC movement a deviation is added to its own region's preference,
\\h_r + \epsilon_r\\, before the gradient \\d = h_i - h_j\\ is taken
along each edge. One deviation therefore moves every edge that touches
that region, raising the rates into it and lowering the rates out of it,
and the perturbed taxis is still the gradient of a surface. Only
differences reach the generator, so a constant added to every region's
deviation leaves movement unchanged and the level of the field is held
only by the process error penalty.
