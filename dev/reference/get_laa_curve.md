# Mean length and its spread at a set of real ages

Schnute-form von Bertalanffy growth: `L1` is the mean length at age
`A1`, `L2` the mean length at age `A2`, and `K` the growth rate, so
\$\$L\_\infty = L_1 + \frac{L_2 - L_1}{1 - e^{-K(A_2 - A_1)}}\$\$ and
\\L(x) = L\_\infty + (L_1 - L\_\infty) e^{-K(x - A_1)}\\ for real ages
\\x \ge A_1\\. Below `A1` growth is linear from `L0` at age zero, \\L(x)
= L_0 + (x / A_1)(L_1 - L_0)\\, the linear phase. `L2_asymptote = 1`
reads `L2` as the asymptote directly, with no second reference age to
solve it from.

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

  Numeric vector of real ages (data, not parameters).

- L0:

  Length at age zero, the anchor of the linear phase.

- L1, L2, K, CV1, CV2:

  Growth parameters, natural scale, possibly AD.

- A1, A2:

  Reference ages for `L1` and `L2`. Ignored for the asymptote under
  `L2_asymptote`, though `A2` still bounds the CV interpolation.

- cv_type:

  Integer, 0 interpolate the CV on length, 1 scale by age.

- sd_type:

  Integer, 0 the CV parameters scale the mean, 1 they are SDs.

- A2_cv:

  Age at and above which `CV2` applies. Defaults to `A2`. Under
  `L2_asymptote` there is no second reference age, so
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md)
  sets `A2` to the accumulator age and the interpolation runs to there.

- rho:

  Richards coefficient, natural scale, possibly AD. One (the default) is
  the von Bertalanffy curve.

- cv_ref:

  Optional vector of the coefficient of variation at each element of
  `x`, used in place of the one this curve implies. Holds the spread at
  age at a reference year's while the mean moves, which is the
  convention for a time-varying growth curve.

- L2_asymptote:

  Integer, 0 (default) solves \\L\_\infty\\ from `L1` and `L2` at their
  reference ages, 1 reads `L2` as \\L\_\infty\\ itself. Set from
  `growth_A2 = "Linf"` in
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

## Value

List with `L` (mean length), `sd` (spread), `Linf` and `cv`.

## Details

With a Richards coefficient `rho` other than one the curve is the
Richards generalization, which applies the same form to the lengths
raised to that power, \$\$L(x)^\rho = L\_\infty^\rho + (L_1^\rho -
L\_\infty^\rho) e^{-K(x - A_1)}\$\$ with \\L\_\infty^\rho = L_1^\rho +
(L_2^\rho - L_1^\rho) / (1 - e^{-K(A_2 - A_1)})\\ when `A2` is a real
age. `rho = 1` is the von Bertalanffy curve.

The coefficient of variation is `CV1` below `A1`, `CV2` at and above
`A2`, and in between interpolates linearly on mean length
(`cv_type = 0`) or on age (`cv_type = 1`). The spread is `CV * L` under
`sd_type = 0` and the parameter itself under `sd_type = 1`.
