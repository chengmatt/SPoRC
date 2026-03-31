# Do Population Projections

Do Population Projections

## Usage

``` r
Do_Population_Projection(
  n_proj_yrs = 2,
  n_regions,
  n_ages,
  n_sexes,
  sexratio,
  n_fish_fleets,
  do_recruits_move = 0,
  recruitment,
  terminal_NAA,
  terminal_NAA0,
  terminal_F,
  natmort,
  WAA,
  WAA_fish,
  MatAA,
  fish_sel,
  Movement,
  f_ref_pt = NULL,
  b_ref_pt = NULL,
  HCR_function = NULL,
  recruitment_opt = "inv_gauss",
  fmort_opt = "HCR",
  t_spawn,
  bh_rec_opt = NULL,
  init_F = 0
)
```

## Arguments

- n_proj_yrs:

  Number of projection years

- n_regions:

  Number of regions

- n_ages:

  Number of ages

- n_sexes:

  Number of sexes

- sexratio:

  Array of recruitment sex ratio (n_regions, n_proj_yrs, n_sexes)

- n_fish_fleets:

  Number of fishery fleets

- do_recruits_move:

  Whether recruits move (0 == don't move, 1 == move)

- recruitment:

  Recruitment matrix dimensioned by n_regions, and n_yrs that we want to
  summarize across, or condition our projection on

- terminal_NAA:

  Terminal Numbers at Age dimensioned by n_regions, n_ages, n_sexes

- terminal_NAA0:

  Terminal Unfished Numbers at Age dimensioned by n_regions, n_ages,
  n_sexes

- terminal_F:

  Terminal fishing mortality rate, dimensioned by n_regions,
  n_fish_fleets

- natmort:

  Natural mortality, dimensioned by n_regions, n_proj_yrs, n_ages,
  n_sexes

- WAA:

  Weight at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes

- WAA_fish:

  Weight at age for the fishery, dimensioned by n_regions, n_proj_yrs,
  n_ages, n_sexes, n_fish_fleets

- MatAA:

  Maturity at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes

- fish_sel:

  Fishery selectivity, dimensioned by n_regions, n_proj_yrs, n_ages,
  n_sexes, n_fish_fleets

- Movement:

  Movement, dimensioned by n_regions, n_regions, n_proj_yrs, n_ages,
  n_sexes

- f_ref_pt:

  Fishing mortality reference point dimensioned by n_regions and
  n_proj_yrs

- b_ref_pt:

  Biological reference point dimensioned by n_regions and n_proj_yrs

- HCR_function:

  Function describing a harvest control rule. The function should always
  have the following arguments: x, which represents SSB, frp, which
  takes inputs of fishery reference points, and brp, which takes inputs
  of biological reference points. Any additional arguments should be
  specified with defaults or hard coded / fixed within the function.

- recruitment_opt:

  Recruitment simulation option, where options are "inv_gauss", which
  simulates future recruitment based on the the recruitment values
  supplied using an inverse gaussian distribution, "mean_rec", which
  takes the mean of the recruitment values supplied for a given region,
  and "zero", which assumes that future recruitment does not occur

- fmort_opt:

  Fishing mortality option. Choices are:

  \* \*\*"HCR"\*\* – Applies the user-supplied \`HCR_function\` using
  region-specific SSB, F reference point, and biomass reference point.

  \* \*\*"HCR_global"\*\* – Applies the \`HCR_function\` using global
  SSB (summed across regions) and a global biomass reference point (sum
  of the region-specific biomass reference points). Each region's
  biomass reference point should be defined individually; the function
  performs the summation.

  \* \*\*"Input"\*\* – Uses user-supplied projected fishing mortality
  values directly.

- t_spawn:

  Fraction time of spawning used to compute projected SSB

- bh_rec_opt:

  A list object containing the following arguments:

  recruitment_dd

  :   A value (0 or 1) indicating global (1) or local density dependence
      (0). In the case of a single region model, either local or global
      will give the same results

  rec_lag

  :   A value indicating the number of years lagged that a given year's
      SSB produces recruits

  R0

  :   The virgin recruitment parameter

  Rec_Prop

  :   Recruitment apportionment values. In a single region model, this
      should be set at a value of 1. Dimensioned by n_regions

  h

  :   Steepness values for the stock recruitment curve. Dimensioned by
      n_regions

  WAA

  :   A weight-at-age array dimensioned by n_regions, n_ages, and
      n_sexes, where the reference year should utilize values from the
      first year

  MatAA

  :   A maturity at age array dimensioned by n_regions, n_ages, and
      n_sexes, where the reference year should utilize values from the
      first year

  natmort

  :   A natural mortality at age array dimensioned by n_regions, n_ages,
      and n_sexes, where the reference year should utilize values from
      the first year

  SSB

  :   All SSB values estimated from a given model, dimensioned by
      n_regions and n_yrs

- init_F:

  Scalar of initial F value to apply for deriving beverton holt
  recruitment; default is set at 0.

## Value

A list containing projected F, catch, SSB (and dynamic unfished), and
Numbers at Age (and dynamic unfished). (Objects are generally
dimensioned in the following order: n_regions, n_yrs, n_ages, n_sexes,
n_fleets)

## See also

Other Reference Points and Projections:
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/reference/Get_Reference_Points.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/reference/get_key_quants.md)
