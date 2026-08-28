# Name the split margins of one at-age stream

An at-age observation is stored over regions and sexes whether or not it
is reported that way, and the Type code names which of those margins the
fleet reports separately. The vocabulary is the composition one, so a
model stating both kinds of data states them the same way.

## Usage

``` r
do_at_age_type_setup(
  input_list,
  type,
  stream,
  fleet_field,
  use_field,
  pop = FALSE
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- type:

  Character, one of `"agg"`, `"spltRaggS"`, `"aggRspltS"` or
  `"spltRspltS"`, either one setting for every fleet or one per fleet.

- stream:

  Stream tag naming the data element, e.g. `"CatchAA"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- use_field:

  Name of the use array for this stream.

- pop:

  Logical. `TRUE` for the population-specific stream.

## Value

`input_list` with `$data$<stream>_Type` set.

## Details

A margin the fleet sums over carries its observation in slot one, and a
use flag anywhere else on that margin is refused rather than quietly
ignored.
