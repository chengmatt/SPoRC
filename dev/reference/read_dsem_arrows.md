# Read dsem arrow and lag notation

One arrow per line. `"from -> to, lag, name, start"` is a path: series
`from` in year `t - lag` affects series `to` in year `t`. A lag may be
negative, which reads `from` in a later year rather than an earlier one.
On a time axis that says the future affects the present, so it is for
axes that are not time, such as ages or length bins.
`"a <-> a, 0, name, start"` is the innovation sd of series `a`, and
`"a <-> b, 0, name, start"` a covariance term between two innovations. A
name of `NA` fixes the arrow at `start`, and a name used on several
arrows is one shared parameter.

## Usage

``` r
read_dsem_arrows(
  dsem_arrows,
  variables,
  covs = NULL,
  mod_var_logscale = FALSE,
  variance = "conditional"
)
```

## Arguments

- dsem_arrows:

  Arrow lines, one per element of a character vector or one per line of
  a string. `#` starts a comment.

- variables:

  Character vector of series names the arrows may use.

- covs:

  Character vector of series groups whose innovations are allowed to
  covary, each group written as one string, for example `"a, b"`. A
  covariance term is added for every pair in a group. Series with no sd
  line of their own get one added.

- mod_var_logscale:

  Whether a moderated sd or covariance is the exponential of its series.
  `FALSE` (default) reads it on the natural scale.

- variance:

  What an sd line means. `"conditional"` (default) makes it the
  innovation sd, so a series' spread is that plus whatever its paths
  add: \\s / \sqrt{1 - \rho^2}\\ under a self path \\\rho\\, which
  increases without bound under a random walk. `"diagonal"` and
  `"marginal"` make it the marginal sd of the series, with the
  innovation sd solved for year by year
  ([`get_dsem_matrices`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_matrices.md)).
  The two are the same without covariance lines. With them, `"diagonal"`
  solves as if the innovations were independent, so a cell lands near
  its sd line, and `"marginal"` keeps the innovation correlations the
  lines imply and lands every cell on it exactly. These are `dsem`'s
  `constant_variance` settings. Both forms need stationary paths and do
  not allow a series with an sd of zero or a moderated sd.

## Value

List with `arrows` (one row per arrow: type, from, to, lag, name, start,
par, from_idx, to_idx, mod_idx), `variables`, `beta_names` (paths and
covariances, natural scale), `ln_sd_names` (sds, log scale),
`mod_var_logscale` and `series_order`, the order the series can be drawn
in within a year (`NULL` when same-year paths form a loop), and
`project_k`, the series whose sd is fixed at zero, worked out from what
points into them, and `variance`, what the sd lines mean.
