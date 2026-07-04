# Map fishery selectivity process error hyperparameters

Constructs the factor map for `fishsel_pe_pars`, which contains the
variance and correlation hyperparameters governing continuous
time-varying selectivity. The set of active parameters depends on the
time-variation type (`cont_tv_fish_sel`): iid/random-walk forms use up
to 2 parameters (log-sigma); 3D GMRF forms use up to 4 (partial
correlations for age, year, cohort dimensions plus log-sigma); the 2D
AR1 form uses 3 (bin AR1, year AR1, log-sigma). Correlation components
can be selectively suppressed via `corr_opt_semipar`.

## Usage

``` r
do_fishsel_pe_pars_mapping(
  input_list,
  fishsel_pe_pars_spec,
  corr_opt_semipar,
  bins
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- fishsel_pe_pars_spec:

  Character vector of length `n_fish_fleets`. Options:

  `"est_all"`

  :   Separate hyperparameters per region × sex.

  `"est_shared_r"`

  :   Shared across regions; unique per sex.

  `"est_shared_s"`

  :   Shared across sexes; unique per region.

  `"est_shared_r_s"`

  :   Shared across regions and sexes.

  `"est_shared_f_x"`

  :   Copy hyperparameters from fleet `x`.

  `"fix"` or `"none"`

  :   All parameters fixed (mapped to `NA`).

- corr_opt_semipar:

  Character vector of length `n_fish_fleets` specifying which
  correlation components to suppress for semi-parametric models. Valid
  values per fleet: `NA` (no suppression), `"corr_zero_y"`,
  `"corr_zero_b"`, `"corr_zero_y_b"`, `"corr_zero_c"`,
  `"corr_zero_y_c"`, `"corr_zero_b_c"`, `"corr_zero_y_b_c"`. Cohort
  options (`"corr_zero_c"`, etc.) are only valid for 3D GMRF forms and
  will error if applied to the 2D AR1 (`cont_tv_fish_sel == 5`).

- bins:

  Number of selectivity bins

## Value

The input `input_list` with `$map$fishsel_pe_pars` set to a factor
vector. Index numbering is reset after any correlation suppression to
maintain contiguous integer indices.

## Details

Fleet sharing (`"est_shared_f_x"`) is handled in a second pass.
