# Map catch-at-age observation error from a key matrix

The key is an integer matrix `[n_ages, n_fish_fleets]` in which equal
entries share a parameter and `NA` excludes one. This is the key matrix
convention ICES age-structured assessments use for coupling. One
structure covers every sharing pattern that would otherwise need its own
spec string: one standard deviation per age, standard deviations by age
group, or one for the whole fleet.

## Usage

``` r
do_sigmaCAA_mapping(input_list, key, spec, starting_values = list())
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- key:

  Integer matrix `[n_ages, n_fish_fleets]`, or `NULL` for one parameter
  per fleet shared across ages. A season needing its own observation
  error is a separate fleet, the same way anything else needing its own
  selectivity or catchability is.

- spec:

  `"est"` or `"fix"`.

- starting_values:

  Named list from the caller's `...`, read for `ln_sigmaCAA`.

## Value

`input_list` with `$par$ln_sigmaCAA` and `$map$ln_sigmaCAA` set.

## Details

Coupled parameters are checked against the observations that inform
them. A standard deviation with a single observation is not merely
poorly determined: the likelihood is unbounded, since the standard
deviation collapses onto whatever residual the model can fit exactly and
the `log(sigma)` term runs to negative infinity. The optimiser reports
convergence either way.
