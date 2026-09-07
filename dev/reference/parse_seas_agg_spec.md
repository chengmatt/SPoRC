# Parse a per-fleet seasonal aggregation specification

Seasonal models can report an observation once a year rather than once a
season. `"spltSeas"` fits the observation in the season it sits in,
which is what every data source did before this setting existed.
`"aggSeas"` sums the model prediction over every season of the year and
fits it against a single observation.

## Usage

``` r
parse_seas_agg_spec(spec, arg_name, n_fleets)
```

## Arguments

- spec:

  Character vector of length `n_fleets`, or a single value given to
  every fleet. The resolved codes `0` and `1` are also taken, so a
  fitted model's settings can be handed straight back to an operating
  model. `NULL` leaves every fleet at `"spltSeas"`.

- arg_name:

  Name of the argument being parsed, used in error messages.

- n_fleets:

  Number of fleets the vector must cover.

## Value

Integer vector of length `n_fleets`, `0` for `"spltSeas"` and `1` for
`"aggSeas"`.

## Details

Under `"aggSeas"` the observation still lives in the season it was
placed in, so exactly one season per region and year may be turned on in
the matching `Use` array. That season is where the likelihood, the
residual and the reported negative log likelihood all land; the
prediction it is compared against is the whole year.
