# Map selectivity deviation parameters (fishery, retention, or survey)

Constructs the factor map for continuous time-varying selectivity
deviations (`ln_fishsel_devs`, `ln_retsel_devs`, or `ln_srvsel_devs`)
across region, year, bin, sex, and fleet. For iid/random-walk forms,
active bins are governed by the fitted selectivity model's parameter
count; for 3D GMRF/2D AR1 forms, every age bin is active. Bin groupings
(`sel_devs_shared_bins`, used by `"est_shared_b"` and its combinations)
apply wherever the deviations are indexed by bin: the GMRF and AR1
forms, and the non-parametric selectivity forms under iid or a random
walk. Fleet sharing (`"est_shared_f_x"`) is handled in a second pass.

## Usage

``` r
do_sel_devs_mapping(
  input_list,
  sel_devs_spec,
  sel_devs_shared_bins,
  bins,
  dont_est_dev_first = NULL,
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

- dont_est_dev_first:

  Integer vector `[n_<fleet_field>]` of 0/1, or `NULL`. Where `1`, the
  deviations of year one are dropped from the map so the walk starts in
  year two and the fixed selectivity parameters hold year one. Refused
  for the GMRF and 2D AR1 forms, whose deviations are a field rather
  than a walk anchored at year one.

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

  Character. Stub for the usage-indicator fields read by
  [`sel_has_data`](https://chengmatt.github.io/SPoRC/dev/reference/sel_has_data.md):
  `Use<use_field>` and its at-age counterpart `Use<use_field>AA`, each
  with a `_pop` variant. `"Catch"` for `"fish"`/`"ret"`; `"SrvIdx"` for
  `"srv"`.

- fleet_label:

  Character. Used only in the collected setup message.

## Value

The input `input_list` with `$map$ln_<prefix>sel_devs` set to a factor
vector, and `$data$map_ln_<prefix>sel_devs` set to the equivalent
integer array.

## Details

Serves fishery, retention, and survey selectivity, selected by `prefix`
exactly as in
[`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md).
