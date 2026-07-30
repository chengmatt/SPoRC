# Do Population Projections

Projects population dynamics forward in time under alternative
recruitment and fishing mortality scenarios. The model initializes from
terminal assessment quantities and advances numbers-at-age through
recruitment, seasonal movement, mortality, ageing, and harvest control
rules across multiple seasons and years.

## Usage

``` r
Do_Population_Projection(
  n_proj_yrs = 2,
  n_pop,
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
  dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
  natmort,
  natal_region,
  WAA,
  WAA_fish,
  MatAA,
  fish_sel,
  ret_sel = array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes,
    n_fish_fleets)),
  Movement,
  sgl_seas_spawning_movement,
  stray_rate,
  f_ref_pt = NULL,
  b_ref_pt = NULL,
  HCR_function = NULL,
  recruitment_opt = "inv_gauss",
  fmort_opt = "HCR",
  t_spawn,
  bh_rec_opt = NULL,
  n_seas = 1,
  seasdur = rep(1/n_seas, n_seas),
  spawn_seas = 1,
  rec_seas_prop = {
     rec_seas_prop = array(0, dim = c(n_pop, n_seas))
    
rec_seas_prop[] <- 1/n_seas
     rec_seas_prop
 },
  Mrate = NULL,
  move_timing = 0
)
```

## Arguments

- n_proj_yrs:

  Integer. Number of projection years.

- n_pop:

  Integer. Number of populations (may exceed regions when natal homing
  is modeled).

- n_regions:

  Integer. Number of spatial regions.

- n_ages:

  Integer. Number of age classes including the plus group.

- n_sexes:

  Integer. Number of sexes.

- sexratio:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_sexes\]\`. Recruitment sex
  ratio used to allocate projected recruits by sex.

- n_fish_fleets:

  Integer. Number of fishing fleets.

- do_recruits_move:

  Integer (0 or 1). Whether age-1 recruits are subject to movement.
  Default = 0.

- recruitment:

  Array \`\[n_pop, n_regions, n_yrs\]\`. Historical recruitment used to
  condition stochastic projection options.

- terminal_NAA:

  Array \`\[n_pop, n_regions, n_seas, n_ages, n_sexes\]\`. Fished
  numbers-at-age in the terminal assessment year.

- terminal_NAA0:

  Array \`\[n_pop, n_regions, n_seas, n_ages, n_sexes\]\`. Unfished
  numbers-at-age in the terminal assessment year.

- terminal_F:

  Array \`\[n_regions, n_seas, n_fish_fleets\]\`. Terminal fishing
  mortality; sets F in projection year 1 and defines the seasonal F
  ratios applied in subsequent years.

- dmr:

  Array `[n_regions, n_seas, n_fish_fleets]`. Discard mortality rate.
  Default behavior is no discard mortality (`dmr = 0`). When combined
  with `ret_sel = 1`, this implies no discarding within a given fleet
  (all catch is retained).

- natmort:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_ages, n_sexes\]\`. Annual
  natural mortality-at-age, scaled internally by season duration.

- natal_region:

  Integer vector \`\[n_pop\]\`. Natal region for each population.

- WAA:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes\]\`.
  Weight-at-age used in spawning biomass calculations.

- WAA_fish:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes,
  n_fish_fleets\]\`. Fishery weight-at-age used in catch biomass
  calculations.

- MatAA:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes\]\`.
  Maturity-at-age.

- fish_sel:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes,
  n_fish_fleets\]\`. Fishery selectivity-at-age.

- ret_sel:

  Array
  `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]`.
  Retention selectivity-at-age. Default behavior corresponds to full
  retention (`ret_sel = 1`), meaning all captured fish are retained
  unless otherwise specified.

- Movement:

  Array \`\[n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages,
  n_sexes\]\`. Seasonal movement transition matrices.

- sgl_seas_spawning_movement:

  Array \`\[n_pop, n_regions, n_regions, n_proj_yrs, n_ages,
  n_sexes\]\`. Spawning movement matrix applied when \`n_seas = 1\` and
  \`n_pop \> 1\` to redistribute fish to natal grounds prior to SSB
  calculation.

- stray_rate:

  Array \`\[n_pop, n_proj_yrs\]\`. Per-population stray rate used when
  accumulating effective SSB contributions across populations.

- f_ref_pt:

  Array \`\[n_regions, n_proj_yrs\]\`. Fishing mortality reference point
  (e.g., F_MSY) or fixed input F, depending on \`fmort_opt\`.

- b_ref_pt:

  Array \`\[n_pop, n_regions, n_proj_yrs\]\`. Biomass reference point
  used in harvest control rules.

- HCR_function:

  Function. Harvest control rule with arguments \`x\` (SSB), \`frp\` (F
  reference point), and \`brp\` (B reference point).

- recruitment_opt:

  Character. Recruitment scenario: \`"inv_gauss"\`, \`"mean_rec"\`,
  \`"zero"\`, or \`"bh_rec"\`.

- fmort_opt:

  Character. Fishing mortality scenario: \`"HCR"\`, \`"HCR_global"\`, or
  \`"Input"\`.

- t_spawn:

  Numeric scalar. Fraction of the spawning season elapsed before
  spawning; used for mid-season SSB calculations.

- bh_rec_opt:

  Named list of inputs for deterministic Beverton–Holt recruitment when
  \`recruitment_opt = "bh_rec"\`. This list is passed directly to
  [`Get_Det_Recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Det_Recruitment.md)
  and must contain all required arguments for that function.

  Required elements and their expected dimensions include:

  `R0`

  :   Numeric vector `[n_pop]`. Unfished recruitment.

  `h`

  :   Numeric array `[n_pop, n_regions]`. Steepness.

  `rec_region_prop`

  :   Numeric array `[n_pop, n_regions]`. Recruitment allocation across
      regions (sums to 1 across regions).

  `rec_seas_prop`

  :   Numeric array `[n_pop, n_seas]`. Seasonal recruitment proportions
      (sums to 1 across seasons).

  `SSB`

  :   Numeric array `[n_pop, n_regions, n_yrs]`. Historical spawning
      biomass, to which projected SSB is appended internally.

  `WAA`

  :   Array `[n_pop, n_regions, n_seas, n_ages]`. Weight-at-age.

  `MatAA`

  :   Array `[n_pop, n_regions, n_seas, n_ages]`. Maturity-at-age.

  `natmort`

  :   Array `[n_pop, n_regions, n_ages]`. Natural mortality.

  `Movement`

  :   Array `[n_pop, n_regions, n_regions, n_seas, n_ages]`. Movement
      transition matrices.

  `sgl_seas_spawning_movement`

  :   Array `[n_pop, n_regions, n_regions, n_ages]`. Spawning movement
      (single-season case).

  `stray_rate`

  :   Numeric vector `[n_pop]`. Straying rates.

  `init_F`

  :   Array `[n_regions, n_seas, n_fish_fleets]`. Initial fishing
      mortality.

  `fish_sel`

  :   Array `[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]`. Total
      selectivity.

  `ret_sel`

  :   Array `[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]`.
      Retention selectivity.

  `dmr`

  :   Array `[n_regions, n_seas, n_fish_fleets]`. Discard mortality
      rates.

  `sex_ratio_f`

  :   Numeric array `[n_pop, n_regions]`. Female recruitment proportion.

  Additional scalar inputs include `rec_dd`, `rec_lag`, `n_pop`,
  `n_regions`, `n_ages`, `n_seas`, `spawn_seas`, `seasdur`, `t_spawn`,
  and `do_recruits_move`.

  Spawning biomass used in recruitment is constructed internally by
  combining `bh_rec_opt$SSB` with projected SSB values during the
  simulation.

  `bh_rec_opt$rec_lag = 1` is the classic lagged case: each projection
  year's recruitment is computed up front from the prior year's SSB,
  exactly as `recruitment_opt = "inv_gauss"`/ `"mean_rec"` are.
  `bh_rec_opt$rec_lag = 0` is age-0 recruitment: recruitment for year
  `y` is computed from year `y`'s own SSB once `spawn_seas` is reached
  within that year's season loop, and is inserted no earlier than
  `spawn_seas` (`rec_seas_prop` must be zero for every season before
  `spawn_seas` in that case). Reference points and the seasonal SBPR
  calculation used to get `bh_rec_opt$WAA`/`MatAA`/etc. are unaffected
  by this choice – `rec_lag` only changes which year's SSB feeds the
  Beverton-Holt curve, not the per-recruit math itself.

- n_seas:

  Integer. Number of seasons. Default = 1.

- seasdur:

  Numeric vector \`\[n_seas\]\`. Duration of each season as a fraction
  of the year.

- spawn_seas:

  Integer. Spawning season index.

- rec_seas_prop:

  Array \`\[n_pop, n_seas\]\`. Proportion of annual recruitment entering
  in each season. Must sum to 1 across seasons for each population.

- Mrate:

  Array dimensioned like \`Movement\`, holding the instantaneous
  movement rates (the generator) rather than the realized transition
  fractions. Only read when \`move_timing\` is 1 or 2. \`NULL\`
  (default) is valid for \`move_timing = 0\`, where movement is applied
  as a transition matrix and no generator is needed.

- move_timing:

  Integer. When movement happens relative to mortality within a season.
  \`0\` (default) applies movement first and mortality afterwards; \`1\`
  applies mortality first and movement afterwards; \`2\` runs the two
  continuously and simultaneously, which also switches catch-at-age to
  the spatial Baranov form built on season-integrated abundance. Must
  match the timing used to derive the reference points the projection is
  run against.

## Value

A named list containing projected fishing mortality, catch, spawning
biomass, effective spawning biomass, dynamic unfished biomass, and
numbers-at-age for fished and unfished states.

## Details

Population dynamics are tracked at full resolution over
`[population x region x year x season x age x sex]`. Recruitment is
generated annually and then distributed across seasons using
`rec_seas_prop`, allowing intra-annual timing of recruitment within the
first age class.

Each projection year proceeds as follows when
`recruitment_opt != "bh_rec"` or `bh_rec_opt$rec_lag != 0` (the classic
case):

1.  Annual recruitment is generated and allocated across regions and
    sexes. Seasonal recruitment is then distributed within the first age
    class using `rec_seas_prop`, with additional recruits entering in
    seasons `seas > 1`.

2.  Fishing mortality-at-age is constructed from annual F, seasonal F
    ratios derived from the terminal year, and selectivity.

3.  Movement is applied at each seasonal step via transition matrices.
    Age-1 movement is optional via `do_recruits_move`.

4.  Within-season mortality is applied using exponential decay. At the
    end of the final season, individuals age forward and the plus group
    accumulates survivors.

5.  Spawning biomass is computed in `spawn_seas` using a mid-season
    mortality correction. For natal homing models with a single season,
    spawning movement is applied prior to SSB calculation.

6.  Catch is calculated using the Baranov equation and aggregated to
    biomass using fishery-specific weights.

7.  Fishing mortality for the next year is updated via the specified
    harvest control rule or fixed input.

When `bh_rec_opt$rec_lag == 0` (age-0 recruitment), steps 1 and 5 above
are reordered within `spawn_seas`: movement is applied first, spawning
biomass is computed from the survivor population alone (no new recruits
exist yet), that SSB is used to generate this year's recruitment, and
only then are the recruits inserted (no earlier than `spawn_seas`) -
immediately before mortality/ageing runs for that season, so the new
cohort is carried forward exactly like any other seasonal recruit pulse.
Years `y > 1` generate recruitment this way; year 1 carries the supplied
terminal assessment state forward with no new recruitment event,
matching the classic case.

Effective spawning biomass at each population's natal region aggregates
contributions from all populations, with cross-population contributions
scaled by `stray_rate` and normalised by the number of populations in
each natal region.

When `n_sexes = 1`, spawning biomass is multiplied by 0.5. When
`n_regions = 1`, movement is skipped.

## See also

Other Reference Points and Projections:
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md)
