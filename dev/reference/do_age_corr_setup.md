# Set the correlation structure for one at-age data source

Each data source is configured where its data are configured, so the
catch and discard data sources are set in
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md)
and the index data sources in their own setup functions. The
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

  `"iid"`, `"1dar1"`, `"us"` or `"2dar1"`, either one setting for every
  fleet or one per fleet.

- data_source:

  Data source tag: `"catch"`, `"discard"`, `"fish_idx"` or `"srv_idx"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- use_field:

  Name of the use array for this data source.

- starting_values:

  Named list from the caller's `...`.

- rho_spec:

  Character string controlling how the correlation parameters are
  shared: `"est_all"`, `"fix"`, or `"est_shared_"` followed by any
  combination of `r`, `s` and `f`, gaining `p` for the
  population-specific data sources. `NULL` (the default) takes
  `"est_shared_r_s"`, or `"est_shared_p_r_s"` when `pop`, both of which
  give one correlation per fleet.

- pop:

  Logical. `TRUE` for the population-specific data source.

## Value

`input_list` with the data source's correlation flag and its correlation
parameters set.

## Details

Four structures are available, per fleet. `"iid"` treats ages as
independent. `"1dar1"` correlates them as an AR(1) in age distance, so a
fleet skipping ages is spaced correctly rather than treated as
consecutive. `"us"` estimates an unstructured correlation across ages,
the third structure ICES age-structured assessments offer. `"2dar1"`
correlates over ages and years jointly through a separable AR(1), which
is defined on a complete grid and so requires the fleet's observed ages
and years to form one.

How the correlations are shared follows the package's spec strings
rather than a structure of its own. They sit over region, sex and fleet,
with a leading population dim for the population-specific data sources,
so `"est_shared_r_s"` (the default, `"est_shared_p_r_s"` for the
population form) gives one per fleet, `"est_shared_r_s_f"` a single
value, `"est_all"` a free one per cell, and `"fix"` holds them all. One
spec governs the data source's across-age correlation, its across-year
correlation and its unstructured matrix together, so two fleets sharing
a correlation share a whole matrix under `"us"`. A region, sex or
population a fleet never observes has no parameter whatever the spec
says, which is what holds the unused slots of a summed dim out.
