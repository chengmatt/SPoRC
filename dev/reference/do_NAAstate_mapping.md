# Set up the state-space numbers at age

Builds the `ln_NAA` parameter array and its map, the process error
standard deviations and their blocking index, and the data fields the
dynamics and the penalty read.

## Usage

``` r
do_NAAstate_mapping(
  input_list,
  NAA_re,
  NAA_re_ages,
  NAA_re_years,
  NAA_sigma_spec,
  NAA_re_region = "iid",
  NAA_re_region_spec = "est_all",
  NAA_re_pop = "iid",
  NAA_re_sex = "iid",
  starting_values = list(),
  NAA_sigma_popblk_spec_vals,
  NAA_sigma_regionblk_spec_vals,
  NAA_sigma_yearblk_spec_vals,
  NAA_sigma_ageblk_spec_vals,
  NAA_sigma_sexblk_spec_vals
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- NAA_re:

  Character. `"none"` (default) leaves the numbers at age deterministic.
  `"iid"` gives every active cell an independent Gaussian innovation.

- NAA_re_ages:

  Ages the state is active over, as ages rather than indices. `NULL`
  (default) uses every age from the second onward.

- NAA_re_years:

  Calendar years the state is active over. `NULL` (default) uses every
  year from the second onward.

- NAA_sigma_spec:

  Character, `"est"` or `"fix"`, whether the process error standard
  deviations are estimated.

- NAA_sigma_popblk_spec_vals, NAA_sigma_regionblk_spec_vals,
  NAA_sigma_yearblk_spec_vals, NAA_sigma_ageblk_spec_vals,
  NAA_sigma_sexblk_spec_vals:

  Lists of integer vectors assigning indices to blocks, exactly as the
  `M_*blk_spec_vals` arguments do. Blocking shares the standard
  deviation; it never removes a cell from the state.

## Value

`input_list` with `$par$ln_NAA`, `$par$ln_sigmaNAA`, their maps, and the
data fields `NAA_re`, `n_est_naa_re`, `naa_re_ages`, `naa_re_yrs` and
`naa_sigma_blocks`.

## Details

The state covers ages 2 and older, including the plus group, over years
2 and later. Year 1 at those ages belongs to `ln_InitDevs` and age 1
belongs to `ln_RecDevs` in every year, so the three parameterizations
partition the numbers at age rather than overlapping. The plus group is
inside the state and not optional: it is the only cell whose influence
never decays, so leaving it deterministic gives up the conditional
independence the state is worth having for.

Ages and years must each be a contiguous run. The penalty scores the
active cells as one rectangular slice, and a gap would either leave
scored cells the dynamics never wrote or force a per-cell branch onto
the tape.
