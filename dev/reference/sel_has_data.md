# Is a fleet's selectivity informed by data in some year?

Reads both the aggregated use array and its at-age counterpart. A fleet
fitting catch at age or an index at age carries no aggregated
observations, so keying only off the aggregated array maps its
selectivity off and silently holds it at the starting value.

## Usage

``` r
sel_has_data(data, use_field, r, f)
```

## Arguments

- data:

  Model data list.

- use_field:

  Stub, `"Catch"`, `"Discard"`, `"FishIdx"` or `"SrvIdx"`.

- r, f:

  Region and fleet.

## Value

`TRUE` when any stream for that region and fleet is fit.
