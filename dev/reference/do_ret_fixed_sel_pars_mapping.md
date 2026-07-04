# Map retained fishery selectivity fixed-effect parameters

Constructs the factor map for `ret_fixed_sel_pars` (e.g., \\a\_{50}\\,
\\k\\, \\a\_{max}\\), controlling whether selectivity shape parameters
are estimated independently or shared across regions, sexes, or fleets.
Cells with no catch data (`UseCatch == 0` and `UseCatch_pop == 0`) are
automatically mapped to `NA`.

## Usage

``` r
do_ret_fixed_sel_pars_mapping(
  input_list,
  ret_fixed_sel_pars_spec,
  bins,
  ret_sel_nonpar_est_bins
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- ret_fixed_sel_pars_spec:

  Character vector of length `n_fish_fleets`. Each element specifies the
  estimation structure for one fleet. Options:

  `"est_all"`

  :   Separate parameters per region × sex × block × bin.

  `"est_shared_r"`

  :   Shared across regions; unique per sex × block × bin.

  `"est_shared_s"`

  :   Shared across sexes; unique per region × block × bin.

  `"est_shared_r_s"`

  :   Shared across regions and sexes; unique per block × bin.

  `"est_shared_f_x"`

  :   Copy parameter mapping from fleet `x`. Fleet `x` must not itself
      use `"est_shared_f_y"`.

  `"fix"` / `"fixed_ret_sel_input"`

  :   All parameters fixed (mapped to `NA`).

- bins:

  Number of selectivity bins.

- ret_sel_nonpar_est_bins:

  Optional list specifying bin groupings for non-parametric selectivity.
  Structure is `[[fleet]][[block]]`, where each element is a list of bin
  index vectors defining grouped parameters.

## Value

Updated `input_list` with `$map$ret_fixed_sel_pars` as a factor array.

## Details

Non-parametric selectivity (model `sel_model == 5`) is supported via
`ret_sel_nonpar_est_bins`, which defines bin groupings that are treated
as individual estimated parameters.

Fleet sharing (`"est_shared_f_x"`) is applied in a second pass after all
base mappings are constructed, copying the full parameter index
structure from the reference fleet.
