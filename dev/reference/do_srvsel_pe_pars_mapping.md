# Map survey selectivity process error (variance/correlation) parameters

Internal helper called by
[`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md)
to construct the TMB/RTMB factor map for `srvsel_pe_pars`
`[n_regions × max(max_sel_pars, 4) × n_sexes × n_srv_fleets]`, the
hyperparameters governing the variance and correlation structure of
continuous time-varying survey selectivity. The second dimension is
padded to at least 4 to accommodate the maximum number of process error
parameters across all semi-parametric forms:

- iid or random walk (`cont_tv_srv_sel` 1–2):

  Up to `max_sel_pars` variance parameters (1, 2, or 6 depending on
  selectivity model).

- 3D GMRF (`cont_tv_srv_sel` 3–4):

  4 parameters: age partial correlation (index 1), year partial
  correlation (index 2), cohort partial correlation (index 3), log-sigma
  (index 4).

- 2D AR1 (`cont_tv_srv_sel` 5):

  3 parameters at indices 1, 2, and 4 (bin, year, log-sigma); index 3 is
  always `NA`.

## Usage

``` r
do_srvsel_pe_pars_mapping(
  input_list,
  srvsel_pe_pars_spec,
  corr_opt_semipar,
  bins
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$n_sexes`,
  `$data$cont_tv_srv_sel`, `$data$srv_sel_model`, `$data$UseSrvIdx`,
  `$data$UseSrvIdx_pop`, and `$data$srv_selex_type`.

- srvsel_pe_pars_spec:

  Character vector `[n_srv_fleets]` specifying the sharing structure.
  Same options as
  [`do_srvsel_devs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_srvsel_devs_mapping.md):
  `"est_all"`, `"est_shared_r"`, `"est_shared_s"`, `"est_shared_r_s"`,
  `"est_shared_f_x"`, `"fix"`, `"none"`.

- corr_opt_semipar:

  Character vector `[n_srv_fleets]` or `NULL`. Specifies which
  correlation components to suppress for 3D GMRF or 2D AR1 likelihoods
  by setting the corresponding parameter indices to `NA`. Valid values:
  `NA`, `"corr_zero_y"`, `"corr_zero_b"`, `"corr_zero_y_b"`,
  `"corr_zero_c"` (3D GMRF only), `"corr_zero_y_c"` (3D GMRF only),
  `"corr_zero_b_c"` (3D GMRF only), `"corr_zero_y_b_c"` (3D GMRF only).
  After suppression, remaining non-`NA` parameters are re-indexed
  sequentially.

- bins:

  Number of selectivity bins

## Value

The input `input_list` with `$map$srvsel_pe_pars` set to a factor
vector. Active parameters receive sequential integer indices; fixed,
inactive, and correlation-suppressed parameters are `NA`.

## Details

When `corr_opt_semipar` is non-`NULL`, specified correlation indices are
set to `NA` and the remaining active parameters are re-indexed
sequentially. Fleet-sharing (`"est_shared_f_x"`) is handled in a second
pass. Parameters are automatically fixed for fleets with no continuous
time-variation (`cont_tv_srv_sel = 0`) or no active survey index data.
