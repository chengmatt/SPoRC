# Map an at-age observation error or catchability from a key matrix

The key is an integer matrix `[n_ages, n_fleets]` in which equal entries
share a parameter and `NA` excludes one. This is the key matrix
convention ICES age-structured assessments use for coupling. One
structure covers every sharing pattern that would otherwise need its own
spec string: one parameter per age, one per age group, or one for the
whole fleet.

## Usage

``` r
do_key_mapping(
  input_list,
  key,
  spec,
  par_name,
  fleet_field,
  use_field,
  starting_values = list(),
  pop = FALSE,
  default_shared = TRUE
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- key:

  Integer matrix `[n_ages, n_fleets]`, or `NULL` for the default given
  by `default_shared`. Gains a leading population dimension when `pop`
  is `TRUE`.

- spec:

  `"est"` or `"fix"`.

- par_name:

  Name of the parameter to map, e.g. `"ln_sigmaCAA"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- use_field:

  Name of the use array informing this parameter, e.g. `"UseCatchAA"`.
  An age a fleet never observes is excluded whatever the key says.

- starting_values:

  Named list from the caller's `...`, read for `par_name`.

- pop:

  Logical. `TRUE` for the population-specific stream.

- default_shared:

  Logical. When `key` is `NULL`, `TRUE` gives one parameter per fleet
  shared across ages and `FALSE` gives one per age and fleet.
  Catchability defaults to the latter.

## Value

`input_list` with `$par$<par_name>` and `$map$<par_name>` set.

## Details

Shared by every at-age stream: the catch, discard and index observation
errors, the age-specific catchabilities, and their population-specific
counterparts.

Coupled parameters are checked against the observations informing them.
A standard deviation with a single observation is not merely poorly
determined: the likelihood is unbounded, since it collapses onto
whatever residual the model can fit exactly and the `log(sigma)` term
runs to negative infinity. The optimiser reports convergence either way.
