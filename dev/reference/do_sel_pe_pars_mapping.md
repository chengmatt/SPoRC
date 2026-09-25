# Map selectivity process error hyperparameters (fishery, retention, or survey)

Builds the factor map for the variance and correlation hyperparameters
of continuous time-varying selectivity (`fishsel_pe_pars`,
`retsel_pe_pars` or `srvsel_pe_pars`), selected by `prefix` exactly as
in
[`do_fixed_sel_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_fixed_sel_pars_mapping.md).
Which parameters are active follows the time-variation form: the iid and
random walk forms use up to two log sigmas, the 3D GMRF forms up to four
(the partial correlations over age, year and cohort plus a log sigma),
and the 2D AR1 form three (the bin and year correlations plus a log
sigma). Fleet sharing is handled in a second pass.

## Usage

``` r
do_sel_pe_pars_mapping(
  input_list,
  pe_pars_spec,
  corr_opt_semipar,
  bins,
  sel_devs_spec,
  sel_devs_shared_bins,
  prefix,
  fleet_field,
  use_field,
  fleet_label
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- pe_pars_spec:

  Character vector, one entry per fleet: `"est_all"`, `"est_shared_r"`,
  `"est_shared_s"`, `"est_shared_r_s"`, those four with `_b` added,
  which put one standard deviation across every bin the fleet reads,
  `"fix"` or `"none"`, or `"est_shared_f_x"`.

- corr_opt_semipar:

  Character vector, one entry per fleet, of which correlation components
  to suppress under the semi-parametric forms: `NA`, `"corr_zero_y"`,
  `"corr_zero_b"`, `"corr_zero_y_b"`, `"corr_zero_c"`,
  `"corr_zero_y_c"`, `"corr_zero_b_c"` or `"corr_zero_y_b_c"`. The
  cohort options are valid for the 3D GMRF forms only.

- bins:

  Number of selectivity bins.

- sel_devs_spec:

  Character vector, one entry per fleet, the deviation specification
  passed to
  [`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md).
  Read only to recognize which dims the deviations are shared over, and
  `"est_shared_f_x"` resolves to the referenced fleet's.

- sel_devs_shared_bins:

  List of integer vectors grouping the bins that share one estimated
  deviation, as passed to
  [`do_sel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sel_devs_mapping.md).

- prefix:

  `"fish"`, `"ret"` or `"srv"`, which drives the field names
  `cont_tv_<prefix>_sel`, `<prefix>_sel_model`, `<prefix>_selex_type`
  and `<prefix>sel_pe_pars`.

- fleet_field:

  Name of the `$data` field giving the number of fleets:
  `"n_fish_fleets"` for `"fish"` and `"ret"`, `"n_srv_fleets"` for
  `"srv"`.

- use_field:

  Stub for the use indicator fields
  [`sel_has_data`](https://chengmatt.github.io/SPoRC/dev/reference/sel_has_data.md)
  reads, `Use<use_field>` and `Use<use_field>AA` with their `_pop`
  variants: `"Catch"` for `"fish"` and `"ret"`, `"SrvIdx"` for `"srv"`.

- fleet_label:

  Used only in the collected setup message.

## Value

`input_list` with `$map$<prefix>sel_pe_pars` set to a factor vector.

## Details

The hyperparameters have to match the deviation series the likelihood
actually evaluates, which is one per shared group read at the group's
lowest bin and first sex. Under iid or a random walk on a non-parametric
fleet the log sigmas are indexed by bin, so `"est_shared_b"` leaves one
per bin group, and sharing deviations across sexes leaves one set for
the first sex under every form. The rest are fixed.
