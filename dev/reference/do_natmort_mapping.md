# Map natural mortality parameters to a block structure

Constructs the `M_blocks` index array and the `ln_M` factor map used by
the TMB/RTMB objective function to share or fix natural mortality
parameters across population, region, year, season, age, and sex
dimensions. Each unique combination of blocks is assigned a sequential
integer ID; all cells within a block share the same `ln_M` parameter.

## Usage

``` r
do_natmort_mapping(
  input_list,
  M_spec,
  M_popblk_spec_vals,
  M_regionblk_spec_vals,
  M_yearblk_spec_vals,
  M_seasblk_spec_vals,
  M_ageblk_spec_vals,
  M_sexblk_spec_vals
)
```

## Arguments

- input_list:

  Named list containing `$data`, `$par`, and `$map` sublists, as
  constructed by upstream setup functions.

- M_spec:

  Character string controlling whether `ln_M` is estimated or fixed. One
  of:

  `"est_ln_M"`

  :   Freely estimate `ln_M` across all defined blocks.

  `"fix"`

  :   Fix all `ln_M` parameters by mapping them to `NA`.

- M_popblk_spec_vals:

  List of integer vectors assigning population indices to blocks, e.g.,
  `list(1, 2)` for two population-specific blocks or `list(1:2)` for a
  single shared block.

- M_regionblk_spec_vals:

  List of integer vectors assigning region indices to blocks, e.g.,
  `list(1:3, 4:5)` for two region blocks.

- M_yearblk_spec_vals:

  List of integer vectors assigning year indices to blocks, e.g.,
  `list(1:10, 11:30)` for two time periods.

- M_seasblk_spec_vals:

  List of integer vectors assigning season indices to blocks, e.g.,
  `list(1, 2)` for two season-specific rates.

- M_ageblk_spec_vals:

  List of integer vectors assigning age indices to blocks, e.g.,
  `list(1:5, 6:10)` for two age groups.

- M_sexblk_spec_vals:

  List of integer vectors assigning sex indices to blocks. Use
  `list(1:2)` for a sex-invariant block or `list(1, 2)` for sex-specific
  mortality.

## Value

The input `input_list` with two fields updated:

- `$map$ln_M`:

  Factor vector of length equal to `prod(dim(par$ln_M))`. Each element
  is an integer estimation index when `M_spec = "est_ln_M"`, or `NA`
  when `M_spec = "fix"`.

- `$data$M_blocks`:

  Integer array of dimensions
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]` mapping
  each population-region-year-season-age-sex cell to its corresponding
  `ln_M` parameter index.
