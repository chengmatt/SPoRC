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
  natal_region = NULL,
  WAA,
  WAA_fish,
  MatAA,
  fish_sel,
  ret_sel = array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes,
    n_fish_fleets)),
  Movement,
  sgl_seas_spawning_movement = NULL,
  stray_rate = NULL,
  f_ref_pt = NULL,
  b_ref_pt = NULL,
  HCR_function = NULL,
  recruitment_opt = "inv_gauss",
  fmort_opt = "HCR",
  catch_input = NULL,
  catch_fallback_opt = if (fmort_opt == "Catch") "HCR" else fmort_opt,
  catch_terminal_yr = FALSE,
  catch_f_max = 5,
  catch_tol = 1e-06,
  catch_max_iter = 100,
  t_spawn,
  srr_opt = NULL,
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

  Integer vector \`\[n_pop\]\`. Natal region for each population. Only
  read when \`n_pop \> 1\`, so \`NULL\` (the default) is valid for
  single-population models.

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
  calculation. Only read in that case, so \`NULL\` (the default) is
  valid otherwise.

- stray_rate:

  Array \`\[n_pop, n_proj_yrs\]\`. Per-population stray rate used when
  accumulating effective SSB contributions across populations. Only read
  when \`n_pop \> 1\`, so \`NULL\` (the default) is valid otherwise.

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

  Character. Fishing mortality scenario: \`"HCR"\`, \`"HCR_global"\`,
  \`"Input"\`, or \`"Catch"\`. \`"Catch"\` solves each projection year's
  fishing mortality so that realised catch matches \`catch_input\`,
  leaving every other model quantity untouched.

- catch_input:

  Catch targets in biomass, used when \`fmort_opt = "Catch"\` and
  ignored otherwise. Either an array \`\[n_regions, n_proj_yrs\]\` of
  annual targets or \`\[n_regions, n_proj_yrs, n_seas\]\` of seasonal
  ones; which shape is supplied decides what gets solved.
  \`catch_input\[r, y\]\` is the catch removed from region \`r\` during
  projection year \`y\`, indexed the same way \`proj_Catch\` is in the
  returned list rather than carrying the one year lag \`f_ref_pt\` does.
  Targets are totals over populations and fleets, and over seasons too
  in the annual case. A target of 0 sets \`F = 0\` there without a
  solve.

  With annual targets one annual F per region is solved for and split
  over seasons at the terminal year seasonal shares, exactly as the
  other \`fmort_opt\` settings do. With seasonal targets that constraint
  is released and a separate F is solved per region and season; the
  split across fleets within a season still stays at terminal year
  ratios, so fleet specific targets are not supported either way. A
  season the terminal year did not fish has no fleet split to inherit
  and so can take no catch, which is an error rather than a silent zero.

  Not every projection year has to carry a target. Set a year to \`NA\`
  and it falls back to \`catch_fallback_opt\` instead, which is the
  usual shape of catch advice: a year or two of agreed catch followed by
  the harvest control rule. \`NA\` (no target, use the fallback) and
  \`0\` (a target of no fishing) are different things. A year has to be
  all target or all \`NA\` across regions and seasons, since splitting
  one annual fallback F across only some seasons has no defensible
  reading; a partly specified year is an error.

  Column 1 is only used when \`catch_terminal_yr = TRUE\`; see that
  argument.

- catch_fallback_opt:

  Character. Which rule sets F in the projection years \`catch_input\`
  leaves \`NA\`: \`"HCR"\`, \`"HCR_global"\`, or \`"Input"\`. Defaults
  to \`"HCR"\` under \`fmort_opt = "Catch"\` and to \`fmort_opt\` itself
  otherwise, where it is unused. Ignored unless \`fmort_opt = "Catch"\`,
  and its usual inputs (\`f_ref_pt\`, \`b_ref_pt\`, \`HCR_function\`)
  are only needed if some year actually falls back. Mind the indexing
  difference when mixing the two: \`catch_input\[r, y\]\` is the catch
  taken in year \`y\`, but \`f_ref_pt\[r, y\]\` sets F in year \`y +
  1\`, which is the lag the HCR and Input options have always carried.

- catch_terminal_yr:

  Logical. Whether projection year 1, which replays the terminal
  assessment year, is also solved against its catch target rather than
  fished at \`terminal_F\`. Default \`FALSE\`. Set it \`TRUE\` for the
  common assessment case where the terminal year's catch is itself a
  projection because the year is not yet complete. Note that this
  overrides the F the assessment estimated for that year, and so changes
  the numbers-at-age entering year 2. Note also that with \`n_seas \>
  1\` the terminal year takes all its seasons from \`terminal_NAA\`
  rather than propagating them, so only the last season's F feeds year
  2; the earlier seasons still take their catch, but do not otherwise
  carry forward.

- catch_f_max:

  Numeric. Upper bound on the F searched when \`fmort_opt = "Catch"\`.
  Default 5. A target that cannot be taken even at this F is
  unreachable, in which case F is capped here, the target is undershot,
  and a warning names the regions involved.

- catch_tol:

  Numeric. Relative catch tolerance for the F solver. Default 1e-6.

- catch_max_iter:

  Integer. Maximum solver iterations per projection year. Default 100.

- t_spawn:

  Numeric scalar. Fraction of the spawning season elapsed before
  spawning; used for mid-season SSB calculations.

- srr_opt:

  Named list of inputs for deterministic stock-recruit recruitment when
  \`recruitment_opt\` is \`"bh_rec"\` or \`"ricker_rec"\`. The curve
  itself is taken from \`recruitment_opt\`, so the same list serves
  both. Formerly \`bh_rec_opt\`.

- bh_rec_opt:

  Deprecated. Former name of \`srr_opt\`; supplying it warns and
  forwards. Supplying both is an error. This list is passed directly to
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
  combining `srr_opt$SSB` with projected SSB values during the
  simulation.

  `srr_opt$rec_lag = 1` is the classic lagged case: each projection
  year's recruitment is computed up front from the prior year's SSB,
  exactly as `recruitment_opt = "inv_gauss"`/ `"mean_rec"` are.
  `srr_opt$rec_lag = 0` is age-0 recruitment: recruitment for year `y`
  is computed from year `y`'s own SSB once `spawn_seas` is reached
  within that year's season loop, and is inserted no earlier than
  `spawn_seas` (`rec_seas_prop` must be zero for every season before
  `spawn_seas` in that case). Reference points and the seasonal SBPR
  calculation used to get `srr_opt$WAA`/`MatAA`/etc. are unaffected by
  this choice, `rec_lag` only changes which year's SSB feeds the
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

A named list of projected quantities. Year index 1 is the terminal
assessment year replayed, so year 2 is the first genuinely projected
year and catch advice for terminal year + 1 is read from index 2.

Several arrays carry a trailing \`n_proj_yrs + 1\` year slot, which is
used inconsistently and is noted per element below. In short:
\`proj_NAA\`, \`proj_NAA0\`, \`proj_F\` and \`proj_F_seas\` fill it, and
\`proj_ZAA\`, \`proj_ret_FAA\` and \`proj_disc_FAA\` leave it at 0.

- `proj_F`:

  Array \`\[n_regions, n_proj_yrs + 1\]\`. Annual fishing mortality by
  region, summed over seasons and fleets. The trailing column holds the
  F the harvest control rule or input would apply in the year after the
  projection; it stays 0 under \`fmort_opt = "Catch"\`, where there is
  no further year to solve a target for.

- `proj_F_seas`:

  Array \`\[n_regions, n_proj_yrs + 1, n_seas\]\`. The same fishing
  mortality broken out by season, so \`rowSums(proj_F_seas\[, y, \])\`
  recovers \`proj_F\[, y\]\` for every \`y\`, including the trailing
  column. This is the only place the answer lives when seasonal catch
  targets are used, since an annual total cannot represent a seasonal
  solve.

- `proj_ret_FAA`:

  Array \`\[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes,
  n_fish_fleets\]\`. Retained fishing mortality-at-age, i.e. the
  component that generates landed catch. Only years \`1:n_proj_yrs\` are
  filled; the trailing year slot stays 0.

- `proj_disc_FAA`:

  Array dimensioned as \`proj_ret_FAA\`, and filled over the same years.
  Discard fishing mortality-at-age, i.e. the component killed but not
  landed, set by \`ret_sel\` and \`dmr\`. Total fishing mortality-at-age
  is the sum of the two.

- `proj_Catch`:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets\]\`.
  Retained catch in biomass. Built from \`proj_ret_FAA\`, so discard
  mortality acts on the population but is not counted here, and this is
  the quantity \`catch_input\` is matched against.

- `proj_SSB`:

  Array \`\[n_pop, n_regions, n_proj_yrs\]\`. Female spawning biomass,
  accumulated in \`spawn_seas\` with the \`t_spawn\` mortality
  correction. Halved when \`n_sexes = 1\`.

- `proj_eff_SSB`:

  Array \`\[n_pop, n_proj_yrs\]\`. Effective spawning biomass at each
  population's natal region, aggregating contributions from every
  population with cross-population terms scaled by \`stray_rate\`. Equal
  to spawning biomass summed across regions when \`n_pop = 1\`.

- `proj_Total_Biom`:

  Array \`\[n_pop, n_regions, n_proj_yrs\]\`. Total biomass over all
  ages and both sexes, accumulated at the same point in the season as
  \`proj_SSB\` and using the same definition the estimation model uses
  for \`Total_Biom\`, so the projected series continues the estimated
  one without a discontinuity at the terminal year.

- `proj_Dynamic_SSB0`:

  Array \`\[n_pop, n_regions, n_proj_yrs\]\`. Spawning biomass the
  population would have carried under the same realised recruitment but
  no fishing, for dynamic depletion.

- `proj_NAA`:

  Array \`\[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages,
  n_sexes\]\`. Fished numbers-at-age held at the start of each season,
  before that season's mortality and ageing. Under the default
  \`move_timing = 0\` movement has already been applied at this point;
  under \`move_timing\` 1 and 2 it has not, since movement is deferred
  into the mortality step. The trailing year slot is filled, and holds
  the numbers carried into the year after the projection ends.

- `proj_NAA0`:

  Array dimensioned as \`proj_NAA\`. The unfished counterpart,
  decremented by natural mortality alone.

- `proj_ZAA`:

  Array \`\[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages,
  n_sexes\]\`. Total mortality-at-age for the season: natural mortality
  scaled by season duration, plus retained and discard fishing mortality
  summed over fleets. Only years \`1:n_proj_yrs\` are filled; the
  trailing year slot stays 0.

- `proj_catch_resid`:

  Array shaped like \`catch_input\`: \`\[n_regions, n_proj_yrs\]\` for
  annual targets, \`\[n_regions, n_proj_yrs, n_seas\]\` for seasonal
  ones. Relative miss on each catch target, \`(realised - target) /
  target\`, and \`NA\` for years carrying no target (including every
  year when \`fmort_opt != "Catch"\`). Should be at or below
  \`catch_tol\` wherever the solve converged, and is worth checking
  directly rather than relying on warnings alone.

## Details

Population dynamics are tracked over
`[population x region x year x season x age x sex]`. Recruitment is
generated annually and then distributed across seasons using
`rec_seas_prop`, allowing intra-annual timing of recruitment within the
first age class.

Each projection year proceeds as follows when
`recruitment_opt != "bh_rec"` or `srr_opt$rec_lag != 0` (the classic
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

When `srr_opt$rec_lag == 0` (age-0 recruitment), steps 1 and 5 above are
reordered within `spawn_seas`: movement is applied first, spawning
biomass is computed from the survivor population alone (no new recruits
exist yet), that SSB is used to generate this year's recruitment, and
only then are the recruits inserted (no earlier than `spawn_seas`) -
immediately before mortality/ageing runs for that season, so the new
cohort is carried forward exactly like any other seasonal recruit pulse.
Years `y > 1` generate recruitment this way; year 1 carries the supplied
terminal assessment state forward with no new recruitment event,
matching the classic case.

Under `fmort_opt = "Catch"` step 7 moves to the front of the following
year instead: the F that lands a catch target depends on that year's own
numbers-at-age, not on the previous year's spawning biomass, so it
cannot be set at the end of the previous year the way the HCR and Input
options are. The year is run repeatedly at trial F values until realised
catch matches the target, then run once more at the accepted F and
committed. Steps 1 to 6 are otherwise unchanged, and no demographic
input is modified: only F moves. Regions are solved jointly rather than
one at a time, because between-season movement (and, under
`move_timing = 2`, the season-integrated abundance) makes each region's
catch depend on the F set in every other region. Seasonal targets are
instead swept forward one season at a time, which is exact because a
season's catch depends only on the F in that season and earlier ones.

Note that the catch solved against is retained catch, matching
`proj_Catch`. Discard mortality still acts on the population through
`dmr` and `ret_sel`, so a fleet that discards will exert more total F
than the target alone implies.

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
