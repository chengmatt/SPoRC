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
  sel_norm_bins = NULL,
  sex_par_offset = NULL,
  sex_scale_offset = NULL,
  sex_apical_offset = NULL,
  sex_scale = NULL,
  nselbins = NULL,
  dbnrml_raw = NULL,
  dbnrml_startbin = NULL
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

- sex_par_offset:

  Integer vector of length `n_fleets` (0/1), or `NULL` for all zero. For
  a fleet flagged 1, the stored fixed-effect selectivity parameters of
  every sex beyond the first are additive offsets on the first sex's
  stored parameters, evaluated as `pars[sex 1] + pars[sex s]` on the
  transformed (usually log) scale. Because most forms exponentiate their
  parameters, an offset of \\\delta\\ makes the sex-\\s\\ natural-scale
  parameter the first sex's value times \\e^{\delta}\\, which is the
  male-offset convention several existing assessments use. An offset
  kept at zero reproduces sex-shared parameters.

- sex_scale_offset:

  Integer vector of length `n_fleets` (0/1), or `NULL` for all zero. For
  a fleet flagged 1, the realized selectivity curve of every sex beyond
  the first is multiplied by
  `exp(sex_scale[region, block, sex, fleet])`, a constant log-scale
  offset on the whole curve. The scaled curve may exceed one, matching
  the convention of a male selectivity offset applied in log space. Not
  meaningful for forms whose scale is standardized away downstream
  (non-parametric forms and the semi-parametric time-varying
  structures); setup refuses those combinations.

- sex_apical_offset:

  Integer vector of length `n_fleets` (0/1), or `NULL`. Where it is 1,
  every sex beyond the first builds its double normal up to
  `exp(sex_scale[region, block, sex, fleet])` rather than to one,
  leaving the first and last bins where their own parameters put them.
  Mutually exclusive with a scale offset on the same fleet.

- sex_scale:

  Array `n_regions x n_blocks x n_sexes x n_fleets` of log-scale curve
  offsets read when `sex_scale_offset` or `sex_apical_offset` flags a
  fleet, or `NULL` when no fleet is flagged.

- nselbins:

  Integer array `n_regions x n_yrs x n_fleets` naming each cell's
  plateau bin (`0` = none), or `NULL`. Bins beyond it are kept at its
  value for the parametric forms (see `n_sel_bins` in
  [`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md));
  the bicubic form's own setup-built plateau is unaffected.

- dbnrml_startbin:

  Integer vector `[n_fleets]` or `NULL`, the bin each fleet's double
  normal anchors its ascending limb at; see
  [`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md).

## Value

List with `sel`, the age-based array
(`n_pop x n_regions x n_yrs_total x n_seas x n_ages x n_sexes x n_fleets`),
and `sel_l`, the length-based array
(`n_regions x n_yrs_total x n_lens x n_sexes x n_fleets`). Only the
array matching `selex_type` is populated; when `selex_type == 1`, `sel`
is filled downstream by multiplying `sel_l` through the size-age
transition matrix.
