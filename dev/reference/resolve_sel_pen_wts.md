# Resolve a modular selectivity penalty weight vector

Converts a user-supplied selectivity penalty weight specification into
the complete named weight vector consumed by
[`Get_sel_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_sel_PE_loglik.md)
and
[`Get_Selex_Smoothness_Penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex_Smoothness_Penalty.md).
Preserves exact backward compatibility with the older single on/off
`cont_tv_*_sel_penalty` flag (which only ever covered `"yr_devs"`,
`"bin_curve"`, `"yr_curve"` – the process-error-deviation-based penalty
terms used when a fleet has a continuous time-varying selectivity model
active) when no explicit weights are supplied. Six additional
`"smooth_*"` terms are included for
[`Get_Selex_Smoothness_Penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex_Smoothness_Penalty.md)'s
modular terms, which operate directly on a fleet's *realized*
selectivity surface and so apply to any selectivity functional form (not
just the bicubic/cubic spline, `Selex_Model == 8`, that originally
motivated them) – for example, a nonparametric (`Selex_Model == 5`)
fleet with discrete time blocks (mirroring an ADMB assessment's
"selectivity change years") can use these same terms to regularize
curvature/stability across those blocks.

## Usage

``` r
resolve_sel_pen_wts(pen_wts, penalty_flag)
```

## Arguments

- pen_wts:

  `NULL`, or a named numeric vector/list giving independent weights for
  any subset of `"yr_devs"`, `"bin_curve"`, `"yr_curve"`,
  `"smooth_bin_curve"`, `"smooth_bin_diff"`, `"smooth_yr_diff"`,
  `"smooth_yr_curve"`, `"smooth_dome"`, `"smooth_mean_center"`. When
  supplied, any term *not* named is set to `0` (disabled) rather than
  falling back to `penalty_flag` – i.e. supplying `pen_wts` opts fully
  into the explicit, modular weight system rather than partially
  retaining the old implicit behavior.

- penalty_flag:

  Logical. Used only when `pen_wts` is `NULL`: the legacy terms
  (`"yr_devs"`, `"bin_curve"`, `"yr_curve"`) are set to `1` if `TRUE`,
  `0` if `FALSE`. The `"smooth_*"` terms are always `0` in this
  fallback, regardless of `penalty_flag`.

## Value

Named numeric vector of length 9:
`c(yr_devs, bin_curve, yr_curve, smooth_bin_curve, smooth_bin_diff, smooth_yr_diff, smooth_yr_curve, smooth_dome, smooth_mean_center)`.
