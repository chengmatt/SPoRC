# Check that a per-fleet setting carries one entry per fleet

These settings are read one fleet at a time inside the mapping loops. A
vector shorter than the fleet count indexes past its end and yields NA,
which the value check downstream then reports as an unrecognized
setting: the caller is told their value is wrong when its length is what
is wrong, and the value named in the message is often one the message
also lists as valid.

## Usage

``` r
check_fleet_spec_length(spec, n_fleets, what, allow_null = FALSE)
```

## Arguments

- spec:

  The setting.

- n_fleets:

  Number of fleets it must cover.

- what:

  Argument name, used in the message.

- allow_null:

  Whether `NULL` is a legal value for this setting. Settings that are
  only read when a feature is switched on pass `TRUE`; settings the
  model always reads leave it `FALSE`, so a missing one is caught here
  rather than further in.

## Value

`invisible(NULL)`. Called for its error.
