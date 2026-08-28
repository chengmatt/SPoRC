# Evaluate one age-disaggregated observation stream

Computes the at-age negative log likelihood for every fleet in one
stream. Observations arrive already transformed by
[`prep_at_age_obs`](https://chengmatt.github.io/SPoRC/dev/reference/prep_at_age_obs.md)
and registered through
[`OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

## Usage

``` r
get_at_age_stream_nLL(
  obs_t,
  use,
  ln_sigma,
  source,
  pop,
  arrays,
  obs_se = NULL,
  sd_form = 0,
  like_type = 0,
  const = 0,
  corr_type = 0,
  trans_rho = 0,
  trans_rho_year = 0,
  us_pars = NULL,
  aa_type = 1
)
```

## Arguments

- obs_t:

  Registered observations for this stream, one element per cell flagged
  in `use`, in [`which()`](https://rdrr.io/r/base/which.html) order, on
  the scale its fleet's likelihood uses.

- use:

  Integer array flagging which cells are fit, dimensioned region by year
  by season by age by sex by fleet, with a leading population dimension
  when `pop` is `TRUE`.

- ln_sigma:

  Log-scale observation error, over age by sex by fleet, with a leading
  population dimension when `pop` is `TRUE`.

- source, arrays:

  Passed to
  [`get_at_age_prediction`](https://chengmatt.github.io/SPoRC/dev/reference/get_at_age_prediction.md).

- pop:

  Logical. `TRUE` for the population-specific stream, whose arrays carry
  a leading population dimension and whose observations are never summed
  over populations.

- obs_se:

  Reported standard errors shaped like `use`, read only by fleets whose
  `sd_form` asks for them.

- sd_form:

  Integer per fleet, see
  [`at_age_obs_sd`](https://chengmatt.github.io/SPoRC/dev/reference/at_age_obs_sd.md).

- like_type:

  Integer per fleet. `0` lognormal, `1` normal.

- const:

  Small constant added inside the log of a lognormal cell, matching the
  aggregated stream's convention.

- corr_type:

  Integer per fleet. `0` `"iid"`, `1` `"1dar1"`, `2` `"us"`, `3`
  `"2dar1"`.

- trans_rho:

  Unconstrained correlation across ages, over sex by fleet.

- trans_rho_year:

  Unconstrained correlation across years, over sex by fleet, read under
  `"2dar1"`.

- us_pars:

  Unconstrained correlation parameters, over pair by sex by fleet, read
  under `"us"`.

- aa_type:

  Integer per fleet naming the split margins, see
  [`at_age_split`](https://chengmatt.github.io/SPoRC/dev/reference/at_age_split.md).

## Value

A list with `nLL` and `pred`, both arrays shaped like `use` and zero
wherever nothing is fit. The predictions are returned so they can be
reported and plotted directly rather than reconstructed.

## Details

Everything that can differ between fleets does: the margins summed over,
the error structure, whether reported standard errors enter, and whether
the density is lognormal or normal. Ages within a cell may be
independent, an AR(1) across ages, or an unstructured correlation
matrix; a fleet may instead correlate over both age and year through a
separable AR(1), which needs the age by year block it is given to be
complete.
