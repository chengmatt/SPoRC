# Validate the type column of a selectivity prior table

Checks the optional `type` column consumed by
[`get_selex_prior`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_prior.md).
A table without the column is all `"par"` rows (the original behavior)
and passes untouched. `"value"` rows are range-checked here because
their `par` column indexes the selectivity grid rather than the
parameter vector, and their `block` must exist in the fleet's block map
to resolve to a year.

## Usage

``` r
validate_selex_prior_types(selex_prior, use_flag, what, sel_blocks, n_bins)
```

## Arguments

- selex_prior:

  Data frame with columns `region`, `fleet`, `block`, `sex`, `par`,
  `mu`, `sd`, and optionally `type`, or `NULL`.

- use_flag:

  Integer (0/1). When `0` the table is returned unchanged and never
  validated, matching how the prior tables are guarded.

- what:

  Character. Name used in error messages.

- sel_blocks:

  Integer array `[region, year, fleet]` mapping model years to
  selectivity blocks.

- n_bins:

  Integer. Size of the grid the data source's selectivity is
  parameterized on (ages or lengths per its selectivity type).

## Value

The validated table, or the input unchanged when `use_flag` is `0`.
