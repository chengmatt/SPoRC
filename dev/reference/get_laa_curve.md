# Mean length and its spread at a set of real ages

Von Bertalanffy growth in Schnute's form, with `L1` the mean length at
age `A1`, `L2` the mean length at age `A2` and `K` the growth rate.
Below `A1` growth is linear from `L0` at age zero, and a Richards
coefficient other than one applies the same form to the lengths raised
to that power. The coefficient of variation is `CV1` below `A1` and
`CV2` at and above `A2`, interpolating between them. Equations are in
the model equations vignette.

## Usage

``` r
get_laa_curve(
  x,
  L0,
  L1,
  L2,
  K,
  CV1,
  CV2,
  A1,
  A2,
  cv_type = 0,
  sd_type = 0,
  A2_cv = NULL,
  rho = 1,
  cv_ref = NULL,
  L2_asymptote = 0
)
```

## Arguments

- x:

  Numeric vector of real ages, data rather than parameters.

- L0:

  Length at age zero, the anchor of the linear phase.

- L1, L2, K, CV1, CV2:

  Growth parameters on the natural scale, possibly AD.

- A1, A2:

  Reference ages for `L1` and `L2`. `A2` is ignored for the asymptote
  under `L2_asymptote`, but still bounds the CV interpolation.

- cv_type:

  Integer, 0 interpolates the CV on length, 1 on age.

- sd_type:

  Integer, 0 has the CV parameters scale the mean, 1 reads them as
  standard deviations.

- A2_cv:

  Age at and above which `CV2` applies, defaulting to `A2`. Under
  `L2_asymptote` there is no second reference age, so
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md)
  sets `A2` to the accumulator age and the interpolation runs to there.

- rho:

  Richards coefficient on the natural scale, possibly AD. One (the
  default) is the von Bertalanffy curve.

- cv_ref:

  Optional vector of the coefficient of variation at each element of
  `x`, used in place of the one this curve implies. It holds the spread
  at age at a reference year's while the mean moves, which is the
  convention for a time-varying growth curve.

- L2_asymptote:

  Integer, 0 (default) solves \\L\_\infty\\ from `L1` and `L2` at their
  reference ages, 1 reads `L2` as \\L\_\infty\\ itself. Set from
  `growth_A2 = "Linf"` in
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

## Value

List with `L` (mean length), `sd` (spread), `Linf` and `cv`.
