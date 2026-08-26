# Evaluate one age-disaggregated observation stream

Walks the observations a stream fits and returns their negative log
likelihood, shaped like the stream's use array so it can be weighted and
reported alongside the aggregated streams.

## Usage

``` r
get_at_age_stream_nLL(
  obs_log,
  use,
  ln_sigma,
  source,
  pop,
  arrays,
  const = 0,
  corr_type = 0,
  rho = 0
)
```

## Arguments

- obs_log:

  Registered log-scale observations for this stream, one element per
  cell flagged in `use`, in
  [`which()`](https://rdrr.io/r/base/which.html) order.

- use:

  Integer array flagging which cells are fit, dimensioned region by year
  by season by age by fleet, with a leading population dimension when
  `pop` is `TRUE`.

- ln_sigma:

  Log-scale observation error, over age and fleet, with a leading
  population dimension when `pop` is `TRUE`.

- source, pop, arrays:

  Passed to
  [`get_at_age_prediction`](https://chengmatt.github.io/SPoRC/dev/reference/get_at_age_prediction.md).

- const:

  Small constant added inside the log, matching the aggregated stream's
  convention.

- corr_type:

  Integer. `0` is `"iid"`, `1` is `"1dar1"`, an AR(1) across ages within
  a cell.

- rho:

  Correlation per fleet, used when `corr_type` is `1`.

## Value

A list with `nLL` and `pred`, both arrays shaped like `use` and zero
wherever nothing is fit. The predictions are returned so they can be
reported and plotted directly rather than reconstructed.

## Details

Observations arrive log-scale and already registered through
[`RTMB::OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).
Registration must happen against the name `getAll` supplied, so the
caller does it: a vector registered under a local name does not link to
the data element, and the objective then diverges from the reported
likelihood.
