# Map initial age-structure deviation parameters

Builds the factor map for `ln_InitDevs`
`[n_pop x n_regions x (n_ages - 1) x n_sexes]`, the log-scale deviations
from the equilibrium initial age structure. Population, region and age
are resolved on a single-sex slice and then expanded across sexes by
`InitDevs_sex_spec`. Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

## Usage

``` r
do_InitDevs_mapping(
  input_list,
  InitDevs_spec,
  rec_dd,
  init_age_devs_shared,
  InitDevs_sex_spec = "est_shared_s"
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`. Requires
  `$data$equil_init_age_strc`, `$data$rec_region_prop_spec`,
  `$data$natal_region` and `$data$rec_dd`.

- InitDevs_spec:

  Sharing structure for `ln_InitDevs`. `"est_shared_pop_r"` gives one
  set of age deviations across every population and region, required
  when `rec_dd = "global"` and `n_regions > 1`. `"est_shared_r"` gives
  one set per population, shared across its regions. `"fix"` holds every
  deviation at zero. `NULL` (default) estimates all independently, which
  is not permitted when `rec_region_prop_spec = 1` and `n_pop > 1`.

- rec_dd:

  Density dependence inherited from
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).
  `"global"` restricts `InitDevs_spec` to `"est_shared_r"` or
  `"est_shared_pop_r"` when `n_regions > 1`.

- init_age_devs_shared:

  Integer vector of length `n_ages - 1` giving the factor level of each
  age position; positions sharing a value share one parameter. Read
  under `equil_init_age_strc = 3`, and respected by
  `InitDevs_spec = "est_shared_r"` (per population, with a population
  offset) and `"est_shared_pop_r"` (globally, no offset).
  `c(1:42, rep(42, 9))` gives 42 free parameters for a 52-age model with
  43 data ages. Default `NULL`.

- InitDevs_sex_spec:

  `"est_shared_s"` (default) maps every sex onto one age curve penalized
  once; `"est_all"` offsets the factor levels per sex so each has its
  own, each penalized. Also builds `data$init_devs_pen_use`, which flags
  one penalized copy of every estimated parameter.

## Value

`input_list` with `$map$ln_InitDevs` set to a factor vector of length
`prod(dim(par$ln_InitDevs))`. Active parameters take sequential
integers; plus-group slots under `equil_init_age_strc = 1` and non-natal
region slots under `rec_region_prop_spec = 1` are `NA`, and their
starting values are reset to `0`.
