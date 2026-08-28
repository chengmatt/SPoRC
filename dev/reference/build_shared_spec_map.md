# Build a factor map from an "est_all"/"fix"/"est_shared\_..." spec string

Convenience wrapper around
[`build_pe_map`](https://chengmatt.github.io/SPoRC/dev/reference/build_pe_map.md)
for the common case of a single spec string (as used by e.g.
`sigmaC_spec`, `sigmaF_spec`) governing a fixed-effect array with no
additional use/fix masking. Validates `spec` against every dimension
combination implied by `dim_abbrev` before building the map, so invalid
specs still fail with the same style of error message as the
hand-written mapping functions.

## Usage

``` r
build_shared_spec_map(
  dims,
  spec,
  dim_abbrev,
  use = NULL,
  what = "parameter",
  min_obs = 2,
  warn_obs = 5
)
```

## Arguments

- dims:

  Named integer vector of array dimensions (see
  [`build_pe_map`](https://chengmatt.github.io/SPoRC/dev/reference/build_pe_map.md)).

- spec:

  Character scalar: `"est_all"`, `"fix"`, or
  `"est_shared_<abbrev>[_<abbrev>...]"`.

- dim_abbrev:

  Named character vector mapping abbreviation to dimension name, given
  in canonical order, e.g.
  `c(r = "region", y = "year", seas = "season", f = "fleet")`. Values
  must match `names(dims)` exactly.

- use:

  Optional array of the same dimensions as `dims`, non-zero where an
  observation informs that cell (e.g. `UseCatch`, `UseSrvIdx`). When
  supplied, the resulting map is checked for observation-error
  parameters that no data can identify. See
  [`check_spec_map_identifiable`](https://chengmatt.github.io/SPoRC/dev/reference/check_spec_map_identifiable.md).
  Defaults to `NULL`, which skips the check and leaves behavior
  unchanged.

- what:

  Character label for the parameter, used in the check's messages.

- min_obs:

  Integer. A group informed by fewer than this many observations raises
  an error. Default `2`.

- warn_obs:

  Integer. A group informed by fewer than this many observations raises
  a warning. Default `5`.

## Value

Factor vector of length `prod(dims)`, suitable for direct assignment to
`input_list$map$<par>`.
