# Map catchability deviations, their sigma and their correlation

A deviation is estimated where its fleet has a continuous catchability
model and that region and fleet hold index observations. Sigma and rho
follow the sharing strings over region and fleet, blanked wherever the
deviations are.

## Usage

``` r
do_q_devs_mapping(
  input_list,
  q_model,
  sigma_q_spec,
  q_rho_spec,
  prefix,
  fleet_field,
  use_field,
  q_re_years = NULL
)
```

## Arguments

- input_list:

  List with `data`, `par` and `map`.

- q_model:

  Integer vector `[n_fleets]` of process error codes.

- sigma_q_spec:

  Sharing string for the deviation sigma, one of `"est_all"`,
  `"est_shared_r"`, `"est_shared_f"`, `"est_shared_r_f"` or `"fix"`.

- q_rho_spec:

  Sharing string for the ar1 correlation, same options.

- prefix:

  `"fish"` or `"srv"`.

- fleet_field:

  Name of the fleet count in `data`.

- use_field:

  Name of the index use array in `data`.

- q_re_years:

  List `[n_fleets]`, each element the model years (as indices into
  `data$years`) that fleet's deviations are estimated over, or `NULL`
  for every year. `NULL` (default) gives every fleet every year.
  Projection years are always estimated. A year left out holds its
  deviation fixed at its starting value (zero unless set by `...`),
  which a random walk or ar1 then steps through like any other fixed
  cell.

## Value

`input_list` with the three maps set.
