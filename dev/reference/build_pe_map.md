# Build a generic sharing/process-error parameter map

Assigns a unique estimation ID to every combination of the "key"
dimensions (those *not* listed in `share_over`), so cells that agree on
the key dimensions get the same ID regardless of their value along the
`share_over` dimensions. This is the general form behind every
`"est_shared_<dims>"` spec used throughout SPoRC's `Setup_*` mapping
functions (e.g. `est_shared_r`, `est_shared_r_seas`, ...): rather than
hand-enumerating one branch per combination of dimensions, callers just
say which dimensions to collapse.

## Usage

``` r
build_pe_map(dims, share_over = character(0))
```

## Arguments

- dims:

  Named integer vector of array dimensions, e.g.
  `c(region = 3, season = 4, fleet = 2)`. Names must be unique and
  non-empty.

- share_over:

  Character vector, subset of `names(dims)`, giving the dimensions
  across which a single parameter is shared. Use `character(0)`
  (default) to estimate a unique parameter per cell (equivalent to
  `"est_all"`), or `names(dims)` to share one value across everything.

## Value

Integer array with `dim = dims` and `names(dim(.)) = names(dims)`,
containing sequential estimation IDs from 1 to the number of unique
groups. No `NA` handling is done here; apply fixing (e.g. `"fix"` or
use-flags) to the result afterwards.
