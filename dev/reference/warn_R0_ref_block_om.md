# Warn when the estimation and operating models start from different R0

The operating model takes `R0_input` by year, and year one of it does
two jobs: it is the equilibrium
[`generate_initial_age_structure`](https://chengmatt.github.io/SPoRC/dev/reference/generate_initial_age_structure.md)
solves from, and it is the R0 that generates year one's recruitment. So
the operating model necessarily starts from the block in force at year
one.

## Usage

``` r
warn_R0_ref_block_om(data, where)
```

## Arguments

- data:

  The model data list.

- where:

  Name of the calling routine, for the message.

## Value

`NULL`, invisibly. Called for the warning.

## Details

The estimation model instead builds its initial age structure from
`R0[R0_ref_block]`. The two agree whenever `R0_ref_block` is the block
covering year one, which is the default and the case that makes physical
sense, since the equilibrium the series starts from is the one the first
block describes. Under any other `R0_ref_block` the estimation model
cannot reproduce the operating model's initial numbers no matter how
well it fits, so a self-test would be measuring that gap rather than the
feature under test.

Only matters under `use_rinit = 0`; with `use_rinit = 1` both sides
initialize from `rinit` and the blocks never enter.
