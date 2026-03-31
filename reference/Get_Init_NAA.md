# Initialize Numbers-at-Age (NAA) for a Population Model

This function generates initial numbers-at-age (NAA) for a structured
population model across regions, sexes, and age classes. It supports
multiple initialization methods, including iterative solution, scalar
geometric series, and matrix geometric series, optionally accounting for
movement and fishing mortality. Initial age deviations can also be
applied.

## Usage

``` r
Get_Init_NAA(
  init_age_strc,
  init_iter,
  n_regions,
  n_sexes,
  n_ages,
  natmort,
  init_F,
  fish_sel,
  R0_r,
  sexratio,
  Movement,
  do_recruits_move,
  NAA,
  ln_InitDevs
)
```

## Arguments

- init_age_strc:

  Integer specifying the initialization method for the age structure: -
  0: Iterative solution to equilibrium - 1: Scalar geometric series
  solution w/o movement in any groups (no movement in all groups) - 2:
  Matrix geometric series solution (generalizes scalar solution with
  movement) - 3: Scalar geometric series solution w/o movement only in
  plus group (no movement in plus groups)

- init_iter:

  Integer; number of iterations to run when \`init_age_strc = 0\`.

- n_regions:

  Integer; number of spatial regions.

- n_sexes:

  Integer; number of sexes.

- n_ages:

  Integer; number of age classes.

- natmort:

  Array of natural mortality rates with dimensions \`\[regions, ages,
  sexes\]\`.

- init_F:

  Numeric; initial fishing mortality applied (0 for unfished
  population).

- fish_sel:

  Array of fishery selectivity with dimensions \`\[regions, ages, sexes,
  fleets\]\`.

- R0_r:

  Numeric vector of recruitment values for each region.

- sexratio:

  Array \`\[regions, sexes\]\` giving the proportion of each sex in
  recruitment.

- Movement:

  Array \`\[regions, regions, ages, sexes\]\` defining movement
  probabilities.

- do_recruits_move:

  Integer; 0 = recruits do not move, 1 = recruits move according to
  \`Movement\`.

- ln_InitDevs:

  Array \`\[regions, ages-1\]\` of log-scale deviations for initial
  numbers-at-age.

## Value

Array of initial numbers-at-age with dimensions \`\[regions, ages,
sexes\]\`.
