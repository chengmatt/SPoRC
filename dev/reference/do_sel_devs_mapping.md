# Map selectivity deviation parameters (fishery, retention, or survey)

Constructs the factor map for continuous time-varying selectivity
deviations (`ln_fishsel_devs`, `ln_retsel_devs`, or `ln_srvsel_devs`)
across region, year, bin, sex, and fleet. For iid/random-walk forms,
active bins are governed by the fitted selectivity model's parameter
count; for 3D GMRF/2D AR1 forms, every age bin is active, optionally
shared via `sel_devs_shared_bins` groupings (`"est_shared_b"` and its
combinations). Fleet sharing (`"est_shared_f_x"`) is handled in a second
pass.

## Usage

``` r
do_sel_devs_mapping(
  input_list,
  sel_devs_spec,
  sel_devs_shared_bins,
  bins,
  prefix,
  fleet_field,
  use_field,
  fleet_label
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- sel_devs_spec:

  Character vector of length `n_<fleet_field>`. Options: `"est_all"`,
  `"est_shared_r"`, `"est_shared_s"`, `"est_shared_r_s"`,
  `"est_shared_b"`, `"est_shared_r_b"`, `"est_shared_b_s"`,
  `"est_shared_r_b_s"`, `"fix"`/`"none"`, or `"est_shared_f_x"`.

- sel_devs_shared_bins:

  List of integer vectors, each defining a group of bins that share a
  single estimated deviation. Required when `sel_devs_spec` includes
  `"est_shared_b"` or its variants.

- bins:

  Number of selectivity bins.

- prefix:

  Character, one of `"fish"`, `"ret"`, or `"srv"`. Drives the
  domain-specific field names: `cont_tv_<prefix>_sel`,
  `<prefix>_sel_model`, `ln_<prefix>sel_devs` (par/map name, no
  underscore before "sel").

- fleet_field:

  Character. Name of the `$data` field giving the number of fleets
  (`"n_fish_fleets"` for `"fish"`/`"ret"`; `"n_srv_fleets"` for
  `"srv"`).

- use_field:

  Character. Stub for the usage-indicator fields: `Use<use_field>` /
  `Use<use_field>_pop`. `"Catch"` for `"fish"`/`"ret"`; `"SrvIdx"` for
  `"srv"`.

- fleet_label:

  Character. Used only in the collected setup message.

## Value

The input `input_list` with `$map$ln_<prefix>sel_devs` set to a factor
vector, and `$data$map_ln_<prefix>sel_devs` set to the equivalent
integer array.

## Details

Consolidates what were previously three near-identical functions,
`do_fishsel_devs_mapping`, `do_retsel_devs_mapping`, and
`do_srvsel_devs_mapping`.
