# Helper function to map natural mortality blocks

This function maps natural mortality (`ln_M`) to a block structure
across region, year, age, and sex dimensions. It assigns unique integer
identifiers to each block defined by the user's specifications, which is
stored in the `M_blocks` array in the input list.

## Usage

``` r
do_M_mapping(
  input_list,
  M_spec,
  M_regionblk_spec_vals,
  M_yearblk_spec_vals,
  M_ageblk_spec_vals,
  M_sexblk_spec_vals
)
```

## Arguments

- input_list:

  A named list object containing model data, parameters, and mapping
  structures.

- M_spec:

  Character string indicating whether to estimate or fix natural
  mortality. Options are:

  - `"est_ln_M"`: Estimate natural mortality parameters.

  - `"fix"`: Fix all natural mortality parameters (requires them to be
    passed in the input list).

- M_regionblk_spec_vals:

  A list of numeric vectors specifying the region indices grouped in
  each block.

- M_yearblk_spec_vals:

  A list of numeric vectors specifying the year indices grouped in each
  block.

- M_ageblk_spec_vals:

  A list of numeric vectors specifying the age indices grouped in each
  block.

- M_sexblk_spec_vals:

  A list of numeric vectors specifying the sex indices grouped in each
  block.

## Value

An updated `input_list` with mapped natural mortality blocks:

- `input_list$map$ln_M` is a factor indicating fixed or estimated
  mortality parameters.

- `input_list$data$M_blocks` is a 4D array (region × year × age × sex)
  with unique block IDs.
