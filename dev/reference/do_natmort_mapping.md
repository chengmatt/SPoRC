# Map natural mortality parameters to a block structure

Builds the `M_blocks` index array and the `ln_M` factor map. Each
combination of blocks takes a sequential integer, and every cell in a
block shares one `ln_M`.

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

  Named list with `$data`, `$par` and `$map`.

- M_spec:

  `"est_ln_M"` estimates `ln_M` across the blocks, `"fix"` maps every
  parameter to `NA`.

- M_popblk_spec_vals:

  List of integer vectors assigning population indices to blocks, e.g.
  `list(1, 2)` or `list(1:2)`.

- M_regionblk_spec_vals:

  List of integer vectors assigning region indices to blocks, e.g.
  `list(1:3, 4:5)`.

- M_yearblk_spec_vals:

  List of integer vectors assigning year indices to blocks, e.g.
  `list(1:10, 11:30)`.

- M_seasblk_spec_vals:

  List of integer vectors assigning season indices to blocks, e.g.
  `list(1, 2)`.

- M_ageblk_spec_vals:

  List of integer vectors assigning age indices to blocks, e.g.
  `list(1:5, 6:10)`.

- M_sexblk_spec_vals:

  List of integer vectors assigning sex indices to blocks, `list(1:2)`
  for one shared rate or `list(1, 2)` for sex-specific mortality.

## Value

`input_list` with `$map$ln_M`, a factor vector of length
`prod(dim(par$ln_M))` holding estimation indices or `NA`, and
`$data$M_blocks`, an integer array
`[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]` giving each
cell's `ln_M` index.
