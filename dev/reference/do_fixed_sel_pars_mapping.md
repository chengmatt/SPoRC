# Map fixed selectivity parameters (fishery, retention, or survey)

Constructs the factor map for fixed-effects selectivity parameters
(`fish_fixed_sel_pars`, `ret_fixed_sel_pars`, or `srv_fixed_sel_pars`)
across region, bin, block, sex, and fleet, handling every parametric
selectivity form (logistic, gamma, exponential, double normal,
asymptotic logistic, bicubic spline) as well as non-parametric
selectivity with user-defined bin groupings. Fleet sharing
(`"est_shared_f_x"`) is handled in a second pass, copying the full
mapping from a reference fleet.

## Usage

``` r
do_fixed_sel_pars_mapping(
  input_list,
  sel_pars_spec,
  bins,
  sel_nonpar_est_bins,
  prefix,
  fleet_field,
  use_field,
  fleet_label
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- sel_pars_spec:

  Character vector of length `n_<fleet_field>`. Options: `"est_all"`,
  `"est_shared_r"`, `"est_shared_s"`, `"est_shared_r_s"`, `"fix"`,
  `"fix_<prefix>_sel_input"`, or `"est_shared_f_x"` (copy from fleet
  `x`).

- bins:

  Number of selectivity bins.

- sel_nonpar_est_bins:

  Optional list `[[fleet]][[block]]` of bin index groupings for
  non-parametric selectivity (model 5).

- prefix:

  Character, one of `"fish"`, `"ret"`, or `"srv"`. Drives the
  domain-specific field names: `use_fixed_<prefix>_sel`,
  `<prefix>_sel_blocks`, `<prefix>_sel_model`,
  `<prefix>_sel_bicubic_binnodes`/`_yrnodes`, `<prefix>_fixed_sel_pars`
  (par/map name).

- fleet_field:

  Character. Name of the `$data` field giving the number of fleets to
  loop over (`"n_fish_fleets"` for both `"fish"` and `"ret"`, since
  retention shares fishery fleets; `"n_srv_fleets"` for `"srv"`).

- use_field:

  Character. Stub for the usage-indicator fields: `Use<use_field>` /
  `Use<use_field>_pop`. `"Catch"` for `"fish"`/`"ret"`; `"SrvIdx"` for
  `"srv"`.

- fleet_label:

  Character. Used only in the collected setup message, e.g.
  `"fishery fleet"` or `"survey fleet"`.

## Value

The input `input_list` with `$map$<prefix>_fixed_sel_pars` set to a
factor vector.

## Details

Consolidates what were previously three near-identical functions,
`do_fish_fixed_sel_pars_mapping`, `do_ret_fixed_sel_pars_mapping`, and
`do_srv_fixed_sel_pars_mapping`.
