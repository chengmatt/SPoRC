# Set up model dimensions

Set up model dimensions

## Usage

``` r
Setup_Mod_Dim(
  years,
  ages,
  lens,
  n_regions,
  n_sexes,
  n_fish_fleets,
  n_srv_fleets,
  n_proj_yrs_devs = 0,
  verbose = FALSE,
  store_config = FALSE
)
```

## Arguments

- years:

  Numeric vector of years.

- ages:

  Numeric vector of age classes.

- lens:

  Numeric vector of length bins; can be set to `1` if length data are
  not modeled.

- n_regions:

  Integer specifying the number of spatial regions.

- n_sexes:

  Integer specifying the number of sexes.

- n_fish_fleets:

  Integer specifying the number of fishery fleets.

- n_srv_fleets:

  Integer specifying the number of survey fleets.

- n_proj_yrs_devs:

  Number of projection years for deviation parameters (ln_RecDevs,
  move_devs, ln_fishsel_devs, ln_srvsel_devs)

- verbose:

  Logical flag indicating whether to print progress messages (default
  `FALSE`).

- store_config:

  Logical flag indicating whether or not to store configuration (default
  `FALSE`).

## Value

A list containing three named elements:

- `data`:

  List of data inputs dimensioned by the model dimensions.

- `parameters`:

  List of model parameters initialized according to dimensions.

- `map`:

  List of parameter mappings for model fitting.

- `config`:

  List of arguments being supplied into the Setup_Mod\_\* functions.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
