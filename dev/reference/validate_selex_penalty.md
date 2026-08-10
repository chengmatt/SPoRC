# Validate a selectivity fixed-effect centering penalty table

Checks the table consumed by
[`get_selex_fixed_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_fixed_penalty.md)
and normalises its `par` column to a list of integer vectors, so a row
may name either one parameter or a whole set.

## Usage

``` r
validate_selex_penalty(selex_penalty, use_flag, what)
```

## Arguments

- selex_penalty:

  Data frame with columns `region`, `fleet`, `block`, `sex`, `par`, and
  `wt`, or `NULL`.

- use_flag:

  Integer (0/1). When `0` the table is returned unchanged and never
  validated, matching how the prior tables are guarded.

- what:

  Character. Name used in error messages.

## Value

The validated table with `par` as a list column, or the input unchanged
when `use_flag` is `0`.
