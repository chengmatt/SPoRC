# Map an estimated index observation error standard deviation

Shared by
[`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md)
and
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
which carry the same parameter under different names and over a
different fleet dimension.

## Usage

``` r
do_sigmaIdx_mapping(input_list, spec, fleet_field, par_name, fleet_map = NULL)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- spec:

  Character scalar. One of `"fix"`, `"est_additive"`, `"est_quadrature"`
  or `"est_replace"`.

- fleet_field:

  Character. Name of the `$data` field giving the number of fleets,
  `"n_fish_fleets"` or `"n_srv_fleets"`.

- par_name:

  Character. Name of the parameter to map, e.g. `"ln_sigmaSrvIdx"`.

- fleet_map:

  Optional integer vector of length `n_fleets`. Fleets sharing a value
  share a parameter, and `NA` fixes that fleet at its starting value.
  Defaults to one free parameter per fleet. Useful when a reference
  assessment estimated some fleets and pinned others at a bound.

## Value

The input `input_list` with `$map$<par_name>` set.

## Details

The parameter is one value per fleet. Index standard errors are reported
per observation, so a year-resolved estimated standard deviation would
place one parameter on one observation and drive the likelihood to
negative infinity; see
[`check_spec_map_identifiable`](https://chengmatt.github.io/SPoRC/dev/reference/check_spec_map_identifiable.md)
for the same problem in the catch and fishing mortality sigmas, which
are dimensioned that way for other reasons.
