# Map survey selectivity fixed-effect parameters

Internal helper called by
[`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md)
to construct the TMB/RTMB factor map for `srv_fixed_sel_pars`
`[n_regions × max_sel_pars × max_sel_blocks × n_sexes × n_srv_fleets]`,
the log-scale fixed-effect parameters of the survey selectivity
functional forms (e.g., \\a\_{50}\\, \\k\\, \\a\_{max}\\). The number of
active parameters per block depends on the selectivity model:
exponential = 1, logistic/gamma = 2, double-normal = 6.

## Usage

``` r
do_srv_fixed_sel_pars_mapping(
  input_list,
  srv_fixed_sel_pars_spec,
  bins,
  srv_sel_nonpar_est_bins
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$n_sexes`,
  `$data$srv_sel_blocks`, `$data$srv_sel_model`, and `$data$UseSrvIdx`.

- srv_fixed_sel_pars_spec:

  Character vector `[n_srv_fleets]` specifying the sharing structure for
  fixed-effect selectivity parameters. One of:

  `"est_all"`

  :   Separate parameters for each region, sex, block, and parameter
      index.

  `"est_shared_r"`

  :   Shared across regions; separate by sex, block, and parameter
      index.

  `"est_shared_s"`

  :   Shared across sexes; separate by region, block, and parameter
      index.

  `"est_shared_r_s"`

  :   Shared across both regions and sexes; separate by block and
      parameter index.

  `"est_shared_f_x"`

  :   Copy the map from fleet `x`. Fleet `x` must not itself use
      `"est_shared_f_x"`.

  `"fix"`

  :   All parameters fixed at starting values (mapped to `NA`).

- bins:

  Number of selectivity bins

## Value

The input `input_list` with `$map$srv_fixed_sel_pars` set to a factor
vector. Active parameters receive sequential integer indices; inactive
region–fleet combinations and fixed parameters are `NA`.

## Details

Parameters are only activated for region–fleet combinations where survey
index data are used (`sum(UseSrvIdx[r,,,f]) > 0`). Fleet-sharing
(`"est_shared_f_x"`) is handled in a second pass after all base mappings
are established, copying the map from the reference fleet.
