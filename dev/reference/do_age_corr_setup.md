# Set the correlation structure for one at-age stream

Each stream is configured where its data are configured, so the catch
and discard streams are set in
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md)
and the index streams in their own setup functions. The
population-specific form carries its own setting rather than borrowing
the aggregated one.

## Usage

``` r
do_age_corr_setup(
  input_list,
  corr,
  stream,
  fleet_field,
  use_field,
  starting_values = list(),
  rho_key = NULL,
  pop = FALSE
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- corr:

  `"iid"`, `"1dar1"`, `"us"` or `"2dar1"`, either one setting for every
  fleet or one per fleet.

- stream:

  Stream tag: `"catch"`, `"discard"`, `"fish_idx"` or `"srv_idx"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- use_field:

  Name of the use array for this stream.

- starting_values:

  Named list from the caller's `...`.

- rho_key:

  Integer matrix `[n_sexes, n_fleets]` coupling the across-age
  correlation, or `NULL` for one per fleet shared across sexes. Equal
  entries share a parameter and `NA` excludes one.

- pop:

  Logical. `TRUE` for the population-specific stream.

## Value

`input_list` with the stream's correlation flag and its correlation
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

Correlations are one per fleet, shared across sexes, unless a key says
otherwise. A sex a fleet never observes carries no parameter.
