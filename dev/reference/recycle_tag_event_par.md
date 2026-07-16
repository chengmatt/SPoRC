# Recycle a per-tag-release-event scalar or vector to the full event count

Internal helper shared by
[`Setup_Sim_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md)
and
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md).
Accepts either a scalar (recycled to all release events) or a vector
already of length `n_events`, and errors on any other length.

## Usage

``` r
recycle_tag_event_par(x, n_events, what)
```

## Arguments

- x:

  Numeric scalar or vector.

- n_events:

  Integer. Number of tag release events (cohorts).

- what:

  Character string used in the error message to identify the offending
  argument.

## Value

Numeric vector of length `n_events`.
