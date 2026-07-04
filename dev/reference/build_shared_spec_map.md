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
build_shared_spec_map(dims, spec, dim_abbrev)
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

## Value

Factor vector of length `prod(dims)`, suitable for direct assignment to
`input_list$map$<par>`.
