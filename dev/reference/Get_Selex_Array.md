# Build the full selectivity array for one selectivity type (retention, fishery, survey)

Evaluates
[`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md)
across regions, years, fleets, and sexes for a single selectivity type,
and mean standardizes the result where the time-varying or
non-parametric form calls for it. Total fishery, retention, and survey
selectivity all share this code path and differ only in which data
arrays are handed in.

## Usage

``` r
Get_Selex_Array(
  selex_type,
  bins,
  sel_blocks,
  sel_model,
  fixed_sel_pars,
  cont_tv_sel,
  ln_seldevs,
  use_fixed_sel,
  sel_input,
  bicubic_Wbin,
  bicubic_Wyr,
  bicubic_binnodes,
  bicubic_yrnodes,
  n_pop,
  n_regions,
  n_yrs,
  n_proj_yrs_devs,
  n_seas,
  n_ages,
  n_lens,
  n_sexes,
  n_fleets,
  bin_devs = NULL,
  bin_dev_bins = NULL,
  sel_norm_bins = NULL
)
```

## Arguments

- selex_type:

  Integer switch (0 = age-based, 1 = length-based).

- bins:

  Numeric vector of bins; `ages` when `selex_type == 0` and `lens` when
  `selex_type == 1`.

- sel_blocks:

  Integer array `n_regions x n_yrs x n_fleets` of selectivity block
  indices.

- sel_model:

  Integer array `n_regions x n_yrs x n_fleets` of selectivity functional
  forms (see
  [`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md)).

- fixed_sel_pars:

  Array `n_regions x n_pars x n_blocks x n_sexes x n_fleets` of
  fixed-effect selectivity parameters.

- cont_tv_sel:

  Integer matrix `n_regions x n_fleets` of continuous time-varying
  forms.

- ln_seldevs:

  Array `n_regions x n_yrs_total x n_bins x n_sexes x n_fleets` of
  selectivity deviations; sliced per fleet internally.

- use_fixed_sel:

  Integer vector of length `n_fleets` (0 = estimate, 1 = use
  `sel_input`).

- sel_input:

  Array
  `n_pop x n_regions x n_yrs x n_seas x n_bins x n_sexes x n_fleets` of
  fixed selectivity values.

- bicubic_Wbin, bicubic_Wyr:

  Bicubic spline bin-node and year-node weight arrays
  (`Selex_Model == 8` only).

- bicubic_binnodes, bicubic_yrnodes:

  Integer arrays `n_regions x n_yrs x n_fleets` giving each block's true
  node counts (`Selex_Model == 8` only).

- n_pop, n_regions, n_yrs, n_proj_yrs_devs, n_seas, n_ages, n_lens,
  n_sexes, n_fleets:

  Model dimensions.

## Value

List with `sel`, the age-based array
(`n_pop x n_regions x n_yrs_total x n_seas x n_ages x n_sexes x n_fleets`),
and `sel_l`, the length-based array
(`n_regions x n_yrs_total x n_lens x n_sexes x n_fleets`). Only the
array matching `selex_type` is populated; when `selex_type == 1`, `sel`
is filled downstream by multiplying `sel_l` through the size-age
transition matrix.
