# Set up the state-space numbers at age

Builds the `ln_NAA` parameter array and its map, the process error
standard deviations and their blocking index, and the data fields the
dynamics and the penalty read. The state covers ages 2 and older
including the plus group, over years 2 and later; year 1 at those ages
belongs to `ln_InitDevs` and age 1 to `ln_RecDevs`. Ages and years must
each be a contiguous run, seasons need not be. A cell is the log numbers
at the start of its season, so season one alone reproduces the annual
state exactly.

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
  NAA_sigma_seasblk_spec_vals,
  NAA_sigma_ageblk_spec_vals,
  NAA_sigma_sexblk_spec_vals,
  NAA_pe_spec = "est_all",
  NAA_re_seasons = "annual",
  NAA_re_season = "iid",
  NAA_re_season_spec = "est_all",
  NAA_re_where = NULL,
  naa_sigma_given = TRUE
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- NAA_re:

  `"none"` (default) leaves the numbers at age deterministic, `"iid"`
  gives every active cell an independent Gaussian innovation.

- NAA_re_ages:

  Ages the state is active over, as ages rather than indices. `NULL`
  (default) uses every age from the second onward.

- NAA_re_years:

  Calendar years the state is active over. `NULL` (default) uses every
  year from the second onward.

- NAA_sigma_spec:

  `"est"` or `"fix"`, whether the process error standard deviations are
  estimated.

- NAA_sigma_popblk_spec_vals, NAA_sigma_regionblk_spec_vals,
  NAA_sigma_yearblk_spec_vals, NAA_sigma_seasblk_spec_vals,
  NAA_sigma_ageblk_spec_vals, NAA_sigma_sexblk_spec_vals:

  Lists of integer vectors assigning indices to blocks, as the
  `M_*blk_spec_vals` arguments. Blocking shares the standard deviation.

- NAA_pe_spec:

  Sharing spec for `NAA_pe_pars`, as documented on
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

- NAA_re_seasons:

  Seasons the state is active over: `"annual"` (default) for season one
  alone, `"all"`, or an integer vector of season indices.

- NAA_re_season:

  `"iid"` (default) or `"us"`, the correlation across the active
  seasons.

- NAA_re_season_spec:

  Sharing spec for the season correlations, taking the same values as
  `NAA_re_region_spec`.

- NAA_re_where:

  Integer matrix `[population, region]`, `1` where the state runs and
  `0` where a population never occupies that region. `NULL` (default)
  gives every cell a state. Cells set to `0` are dropped from the map
  and the penalty, and need the region and population correlations off.

- naa_sigma_given:

  Logical, whether the caller set `NAA_sigma_spec` itself; an explicit
  `"est"` is refused under `NAA_re = "dsem"`.

## Value

`input_list` with `$par$ln_NAA`, `$par$ln_sigmaNAA`, their maps, and the
data fields `NAA_re`, `n_est_naa_re`, `naa_re_ages`, `naa_re_yrs`,
`naa_re_seas` and `naa_sigma_blocks`.
