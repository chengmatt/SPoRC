# Map selectivity process error hyperparameters (fishery, retention, or survey)

Constructs the factor map for the variance/correlation hyperparameters
governing continuous time-varying selectivity (`fishsel_pe_pars`,
`retsel_pe_pars`, or `srvsel_pe_pars`). The set of active parameters
depends on the time-variation type: iid/random-walk forms use up to 2
parameters (log-sigma); 3D GMRF forms use up to 4 (partial correlations
for age, year, cohort dimensions plus log-sigma); the 2D AR1 form uses 3
(bin AR1, year AR1, log-sigma). Correlation components can be
selectively suppressed via `corr_opt_semipar`. Fleet sharing
(`"est_shared_f_x"`) is handled in a second pass.

## Usage

``` r
do_sel_pe_pars_mapping(
  input_list,
  pe_pars_spec,
  corr_opt_semipar,
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

- pe_pars_spec:

  Character vector of length `n_<fleet_field>`. Options: `"est_all"`,
  `"est_shared_r"`, `"est_shared_s"`, `"est_shared_r_s"`,
  `"fix"`/`"none"`, or `"est_shared_f_x"`.

- corr_opt_semipar:

  Character vector of length `n_<fleet_field>` specifying which
  correlation components to suppress for semi-parametric models (`NA`,
  `"corr_zero_y"`, `"corr_zero_b"`, `"corr_zero_y_b"`, `"corr_zero_c"`,
  `"corr_zero_y_c"`, `"corr_zero_b_c"`, `"corr_zero_y_b_c"`). Cohort
  options are only valid for 3D GMRF forms.

- bins:

  Number of selectivity bins.

- prefix:

  Character, one of `"fish"`, `"ret"`, or `"srv"`. Drives the
  domain-specific field names: `cont_tv_<prefix>_sel`,
  `<prefix>_sel_model`, `<prefix>_selex_type`, `<prefix>sel_pe_pars`
  (par/map name, no underscore before "sel").

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

The input `input_list` with `$map$<prefix>sel_pe_pars` set to a factor
vector.

## Details

Serves fishery, retention, and survey selectivity, selected by `prefix`
exactly as in
[`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md).
