# Get Deterministic Recruitment

Computes deterministic recruitment for each region based on either a
mean recruitment model or a Beverton–Holt stock–recruitment
relationship.

## Usage

``` r
Get_Det_Recruitment(
  recruitment_model,
  recruitment_dd,
  y,
  rec_lag,
  R0,
  Rec_Prop,
  h,
  n_regions,
  n_ages,
  WAA,
  MatAA,
  natmort,
  SSB_vals,
  Movement,
  do_recruits_move,
  t_spawn,
  sex_ratio_f,
  init_F,
  fish_sel
)
```

## Arguments

- recruitment_model:

  Integer flag specifying the recruitment model:

  - \`0\` = mean recruitment

  - \`1\` = Beverton–Holt recruitment with steepness

- recruitment_dd:

  Integer flag specifying the scale of density dependence:

  - \`0\` = local (region-specific)

  - \`1\` = global (shared across regions)

- y:

  Current model year (used for SSB lag indexing)

- rec_lag:

  Recruitment lag (number of years between spawning and recruitment)

- R0:

  Virgin or mean recruitment (global scalar)

- Rec_Prop:

  Vector of recruitment proportions by region (used to allocate global
  \`R0\` under local density dependence)

- h:

  Vector of Beverton–Holt steepness values by region

- n_regions:

  Number of spatial regions

- n_ages:

  Number of modeled age classes

- WAA:

  Matrix of weight-at-age by region and age

- MatAA:

  Matrix of maturity-at-age by region and age

- natmort:

  Matrix or vector of natural mortality by region and age

- SSB_vals:

  Matrix of spawning stock biomass (SSB) by region and year

- Movement:

  3D array of movement probabilities between regions by age (\`\[origin,
  destination, age\]\` or alternatively, \`\[n_regions, n_regions,
  age\]\`)

- do_recruits_move:

  Logical or integer flag (0/1) indicating whether recruits move during
  their first year

- t_spawn:

  Fraction of the year at which spawning occurs (used for survival to
  spawning)

- sex_ratio_f:

  Vector of female proportions by region (used to scale initial
  recruits)

- init_F:

  Scalar for initial F value to apply

- fish_sel:

  Array of fishery selectivity of dominant fleet (fleet 1) dimensioned
  by n_regions x n_ages

## Value

A numeric vector of length \`n_regions\` containing deterministic
recruitment values for each region.

## Details

The function returns region-specific deterministic recruitment estimates
based on the chosen recruitment model and density dependence structure.

When \`recruitment_model = 0\`, recruitment is fixed at mean values
(\`R0 \* Rec_Prop\`). When \`recruitment_model = 1\`, Beverton–Holt
recruitment is applied using: \$\$R = \frac{4hR_0SSB}{(1 - h)S_0 + (5h -
1)SSB}\$\$ where \`S_0\` is unfished spawning biomass per recruit,
computed separately for each region (local) or summed across all regions
(global).
