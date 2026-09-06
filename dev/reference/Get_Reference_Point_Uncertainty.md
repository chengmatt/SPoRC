# Confidence intervals for reference points

Attaches uncertainty to
[`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
which returns point estimates only.

## Usage

``` r
Get_Reference_Point_Uncertainty(
  obj,
  SPR_x = NULL,
  t_spawn = 0,
  sex_ratio_f = NULL,
  calc_rec_st_yr = 1,
  rec_age = 1,
  type,
  what,
  n_avg_yrs = 1,
  local_bh_msy_newton_steps = 6,
  is_discard_fleet = NULL,
  sd_rep = NULL,
  method = "delta",
  rel_step = 0.001,
  min_step = 1e-07,
  n_draw = 300,
  seed = NULL,
  level = 0.95,
  extra_quantities = NULL,
  par_subset = NULL
)
```

## Arguments

- obj:

  Fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- SPR_x:

  Numeric. Target spawning potential ratio fraction.

- t_spawn:

  Numeric. Spawning season fraction elapsed before spawning. Default =
  0.

- sex_ratio_f:

  Numeric array `[n_pop, n_regions]`. Defaults to 0.5.

- calc_rec_st_yr:

  Integer. First year of the mean recruitment window. Default = 1.

- rec_age:

  Integer. Recruitment lag in years. Default = 1.

- type:

  Character. `"single_region"` or `"multi_region"`.

- what:

  Character. Reference point method, as in
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md).

- n_avg_yrs:

  Integer. Terminal years averaged for demographic rates. Default = 1.

- local_bh_msy_newton_steps:

  Integer. Newton steps for `"local_MSY"`. Default = 6.

- is_discard_fleet:

  Integer vector `[n_fish_fleets]`. Defaults to all zeros.

- sd_rep:

  Optional `sdreport` to reuse. Default `NULL`.

- method:

  Character. `"delta"`, `"mvn"`, or `"both"`. Default `"delta"`.

- rel_step:

  Numeric. Difference step as a multiple of each parameter's standard
  error. Default 1e-3. Vary it and confirm `d` is unchanged.

- min_step:

  Numeric. Floor on the absolute step. Default 1e-7.

- n_draw:

  Integer. Draws for `"mvn"` and `"both"`. Default 300.

- seed:

  Integer or `NULL`. Seed for the draws. Default `NULL`.

- level:

  Numeric. Confidence level. Default 0.95.

- extra_quantities:

  Optional `function(rep, refpts)` returning a named vector of extra
  positive quantities to propagate through, e.g. stock status.

- par_subset:

  Optional character vector of parameter names to perturb. This is the
  lever for a large model if taking too long. Everything left out is
  asserted to have exactly zero effect, which is a claim about the model
  rather than a shortcut. Make sure to check any subset against
  `method = "mvn"`, whose draws perturb everything.

## Value

A named list:

- `refpts`:

  Delta method estimates, log scale SE, bounds, and CV.

- `mvn`:

  Draw based estimates and quantiles, when requested.

- `d`:

  Sensitivities `[n_quantity, n_par]`.

- `log_cov`:

  Covariance of the log quantities `[n_quantity, n_quantity]`, needed
  for anything combining two of them.

- `draws`:

  Log scale draws `[n_draw, n_quantity]`, when requested.

- `n_par`:

  Parameters perturbed.

- `n_zero_par`:

  How many of those had no effect on any quantity.

- `step`:

  Absolute step used per parameter.

## Details

Write \\\mathbf{p}\\ for everything the fit estimates and \\x^\ast =
\log F^\ast\\ for the reference point. The reference point has no closed
form, so its sensitivity \\d_j = \partial x^\ast / \partial p_j\\ is
taken by central differences of the solved value, using a step of
`rel_step` times each parameter's own standard error. The delta method
then gives \\\mathrm{Var}(x^\ast) = d \Sigma d'\\ with \\\Sigma\\ the
joint covariance of the fit. Intervals are built on the log scale and
exponentiated, so they stay positive and are asymmetric.

`method = "mvn"` skips the linearization and draws parameter vectors
from \\\Sigma\\ instead, re-solving the reference point at each. Use it
to check the delta method when a reference point sits near
\\F\_{crash}\\ or the stock recruit curve is depensatory.

Only estimated quantities contribute. Anything the model holds fixed,
including maturity, fixed \\M\\, fixed steepness, and anything turned
off through `mapping`, gives exactly zero, so an interval can look tight
simply because what drives it was fixed.

## References

Albertsen, C.M. and Trijoulet, V. (2020). Model-based estimates of
reference points in an age-based state-space stock assessment model.
Fisheries Research 230, 105618.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("dusky_rtmb_model")

rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model, SPR_x = 0.4,
                                      type = "single_region", what = "SPR")
rp$refpts

# have stock status through, so the correlation with the reference point is kept
status <- function(rep, refpts) {
  ssb <- rep$SSB[1, 1, dim(rep$SSB)[3]]
  c(SSB_terminal = ssb, status = ssb / as.numeric(refpts$b_ref_pt))
}

rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model, SPR_x = 0.4,
                                      type = "single_region", what = "SPR",
                                      method = "both", extra_quantities = status)
} # }
```
