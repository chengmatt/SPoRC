# Set up mapping for growth

Builds the estimation maps for the growth parameters, for the deviations
of any parameter that varies over time, and for the semi-parametric
surface on mean length at age. Called by
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md)
once the parameter list is populated, since every map is dimensioned off
its parameter.

## Usage

``` r
do_growth_mapping(
  input_list,
  growth_spec,
  growth_fix,
  tv_vals,
  tv_active,
  growth_tv_spec,
  growth_tv_sigma_spec,
  semipar_val,
  growth_semipar_spec,
  semipar_age_idx,
  semipar_yr_idx
)
```

## Arguments

- input_list:

  List containing data, parameter, and map lists, with the growth
  parameters already populated.

- growth_spec:

  Character. How the growth parameters are estimated, one of
  `"est_all"`, `"est_shared_r"`, `"est_shared_s"`, `"est_shared_r_s"` or
  `"fix"`.

- growth_fix:

  Logical vector, one entry per growth parameter, naming any kept at its
  starting value while the others are estimated.

- tv_vals:

  Integer vector, one entry per growth parameter, of the time variation
  each has (0 none, 1 iid, 2 random walk).

- tv_active:

  Matrix `[n_years x n_gpars]` of ones in the years each parameter's
  deviations are estimated in.

- growth_tv_spec:

  Character. How the deviation series are shared across regions and
  sexes.

- growth_tv_sigma_spec:

  Character. `"fix"` holds the process error standard deviations,
  `"est"` estimates them.

- semipar_val:

  Integer code of the semi-parametric form (0 none, 1 iid, 2 random
  walk, 3 `3dmarg`, 4 `3dcond`, 5 `2dar1`).

- growth_semipar_spec:

  Character. Whether the surface's process error parameters are
  estimated.

- semipar_age_idx, semipar_yr_idx:

  Integer vectors of the age and year indices the surface is estimated
  over.

## Value

The input `input_list` with `$map$ln_growth_pars`,
`$map$ln_growth_devs`, `$map$growth_pe_pars` and
`$map$ln_growth_semipar_devs` set, along with the
`$data$map_ln_growth_devs` and `$data$map_ln_growth_semipar_devs`
mirrors the deviation penalties read.
