# Name the split dims of one at-age data source

An at-age observation is stored over regions and sexes whether or not it
is reported that way, and the Type code names which of those dims the
fleet reports separately. The vocabulary is the composition one, so a
model stating both kinds of data states them the same way.

## Usage

``` r
do_at_age_type_setup(
  input_list,
  type,
  data_source,
  fleet_field,
  use_field,
  pop = FALSE
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- type:

  Character. Either bare values, one of `"agg"`, `"spltRaggS"`,
  `"aggRspltS"` or `"spltRspltS"`, given once for every fleet or once
  per fleet; or year and fleet specifications such as
  `"spltRaggS_Year_1-20_Fleet_1"`, in which case every year of every
  fleet needs an entry covering it.

- data_source:

  Data source tag naming the data element, e.g. `"CatchAA"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- use_field:

  Name of the use array for this data source.

- pop:

  Logical. `TRUE` for the population-specific data source.

## Value

`input_list` with `$data$<data source>_Type` set to a numeric matrix
`[n_yrs, n_fleets]` of codes.

## Details

A dim the fleet sums over has its observation in slot one, and a use
flag anywhere else on that dim is refused rather than quietly ignored.

The setting may change part way through a series, given in the same
`Value_Year_x-y_Fleet_f` form the composition data sources take. A bare
value stands for the whole series, either one setting for every fleet or
one per fleet.
