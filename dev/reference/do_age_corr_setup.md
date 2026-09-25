# Set the correlation structure for one at-age data source

Each data source is configured where its data are, so the catch and
discard data sources are set in
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md)
and the index ones in their own setup functions, and the
population-specific form has its own setting rather than borrowing the
aggregated one.

## Usage

``` r
do_age_corr_setup(
  input_list,
  corr,
  data_source,
  fleet_field,
  use_field,
  starting_values = list(),
  rho_spec = NULL,
  pop = FALSE
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- corr:

  Correlation across ages, one setting for every fleet or one per fleet.
  `"iid"` treats ages as independent. `"1dar1"` correlates them as an
  AR(1) in age distance, so a fleet skipping ages is spaced correctly
  rather than treated as consecutive. `"us"` estimates an unstructured
  correlation across ages. `"2dar1"` correlates over ages and years
  jointly through a separable AR(1), which is defined on a complete grid
  and so needs the fleet's observed ages and years to form one.

- data_source:

  Data source tag: `"catch"`, `"discard"`, `"fish_idx"` or `"srv_idx"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- use_field:

  Name of the use array for this data source.

- starting_values:

  Named list from the caller's `...`.

- rho_spec:

  How the correlation parameters are shared: `"est_all"`, `"fix"`, or
  `"est_shared_"` followed by any combination of `r`, `s` and `f`,
  gaining `p` for the population-specific data sources. They sit over
  region, sex and fleet, so `"est_shared_r_s"` gives one per fleet and
  `"est_shared_r_s_f"` a single value. `NULL` (default) takes
  `"est_shared_r_s"`, or `"est_shared_p_r_s"` under `pop`. One spec
  governs the across-age correlation, the across-year correlation and
  the unstructured matrix together, so two fleets sharing under `"us"`
  share a whole matrix. A region, sex or population a fleet never
  observes has no parameter whatever the spec says.

- pop:

  Logical. `TRUE` for the population-specific data source.

## Value

`input_list` with the data source's correlation flag and its correlation
parameters set.
