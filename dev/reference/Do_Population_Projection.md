# Do Population Projections

Projects the population forward under a recruitment and a fishing
mortality scenario, starting from the terminal assessment year and
advancing numbers at age over
`[population x region x year x season x age x sex]` through recruitment,
movement, mortality, ageing and a harvest control rule. Recruitment is
generated annually and spread over seasons by `rec_seas_prop`.

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
  move_timing = 0,
  expm_nsub = 0,
  rec_devs = NULL,
  naa_devs = NULL,
  F_constraint = NULL
)
```

## Arguments

- n_proj_yrs:

  Integer. Number of projection years.

- n_pop:

  Integer. Number of populations, which may exceed the number of regions
  under natal homing.

- n_regions:

  Integer. Number of spatial regions.

- n_ages:

  Integer. Number of age classes including the plus group.

- n_sexes:

  Integer. Number of sexes.

- sexratio:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_sexes\]\` allocating
  projected recruits by sex.

- n_fish_fleets:

  Integer. Number of fishing fleets.

- do_recruits_move:

  Integer (0 or 1). Whether age-1 recruits move. Default 0.

- recruitment:

  Array \`\[n_pop, n_regions, n_yrs\]\` of historical recruitment, used
  to condition the stochastic options.

- terminal_NAA:

  Array \`\[n_pop, n_regions, n_seas, n_ages, n_sexes\]\` of fished
  numbers at age in the terminal assessment year.

- terminal_NAA0:

  As `terminal_NAA`, unfished.

- terminal_F:

  Array \`\[n_regions, n_seas, n_fish_fleets\]\`. Sets F in projection
  year 1 and the seasonal F ratios used in later years.

- dmr:

  Array `[n_regions, n_seas, n_fish_fleets]` of discard mortality rate.
  Default `0`, which with `ret_sel = 1` means a fleet discards nothing.

- natmort:

  Natural mortality, a rate per year in each season, either \`\[n_pop,
  n_regions, n_proj_yrs, n_seas, n_ages, n_sexes\]\` or the same without
  the season dim, which is expanded across seasons. Scaled internally by
  season duration.

- natal_region:

  Integer vector \`\[n_pop\]\` of each population's natal region. Only
  read when \`n_pop \> 1\`, so \`NULL\` (default) is otherwise valid.

- WAA:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes\]\` of
  weight-at-age, used for spawning biomass.

- WAA_fish:

  As `WAA` with a trailing \`n_fish_fleets\` dim, used for catch
  biomass.

- MatAA:

  Array dimensioned like `WAA`, maturity-at-age.

- fish_sel:

  Array \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes,
  n_fish_fleets\]\` of fishery selectivity-at-age.

- ret_sel:

  Array
  `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]`
  of retention selectivity-at-age. Default `1`, full retention.

- Movement:

  Array \`\[n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages,
  n_sexes\]\` of seasonal movement transition matrices.

- sgl_seas_spawning_movement:

  Array \`\[n_pop, n_regions, n_regions, n_proj_yrs, n_ages, n_sexes\]\`
  redistributing fish to natal grounds before SSB. Only read when
  \`n_seas = 1\` and \`n_pop \> 1\`, so \`NULL\` (default) is otherwise
  valid.

- stray_rate:

  Array \`\[n_pop, n_proj_yrs\]\` used when accumulating effective SSB
  across populations. Only read when \`n_pop \> 1\`, so \`NULL\`
  (default) is otherwise valid.

- f_ref_pt:

  Array \`\[n_regions, n_proj_yrs\]\` of the fishing mortality reference
  point or fixed input F, depending on \`fmort_opt\`.

- b_ref_pt:

  Array \`\[n_pop, n_regions, n_proj_yrs\]\` of the biomass reference
  point used in the control rule.

- HCR_function:

  Harvest control rule taking \`x\` (SSB), \`frp\` and \`brp\`. A rule
  that also declares a \`state\` argument (or \`...\`) is given this
  year's population as a named list, holding \`y\`, \`r\`, \`NAA\`,
  \`SSB\`, \`Total_Biom\` and \`Catch\`, so it can be written on more
  than spawning biomass. Rules without it are called as before and no
  state is assembled.

- recruitment_opt:

  Recruitment scenario: \`"inv_gauss"\`, \`"mean_rec"\`, \`"zero"\`, or
  \`"bh_rec"\`.

- fmort_opt:

  Fishing mortality scenario: \`"HCR"\`, \`"HCR_global"\`, \`"Input"\`,
  or \`"Catch"\`, which solves each year's F so realized catch matches
  \`catch_input\` and leaves every other quantity untouched.

- catch_input:

  Catch targets in biomass, read under \`fmort_opt = "Catch"\`. Either
  \`\[n_regions, n_proj_yrs\]\` of annual targets or \`\[n_regions,
  n_proj_yrs, n_seas\]\` of seasonal ones, and the shape decides what is
  solved. \`catch_input\[r, y\]\` is the catch removed in projection
  year \`y\`, indexed as \`proj_Catch\` is rather than with the one year
  lag \`f_ref_pt\` uses. Targets are totals over populations and fleets,
  and over seasons in the annual case; a target of 0 sets \`F = 0\`
  without a solve. Inside the model the F that takes the catch is solved
  numerically first, then refined with two Newton steps the model
  differentiates through, so the F moves with the estimates.

  Annual targets solve one F per region, split over seasons at the
  terminal year shares. Seasonal targets solve an F per region and
  season, with the fleet split within a season still at terminal year
  ratios, so fleet-specific targets are not supported either way. A
  season the terminal year did not fish has no fleet split to inherit
  and can take no catch, which is an error.

  A year set to \`NA\` falls back to \`catch_fallback_opt\`, which is
  the usual shape of catch advice. \`NA\` and \`0\` are different
  things. A year must be all target or all \`NA\` across regions and
  seasons; a partly specified year is an error. Column 1 is read only
  when \`catch_terminal_yr = TRUE\`.

- catch_fallback_opt:

  Which rule sets F in the years \`catch_input\` leaves \`NA\`:
  \`"HCR"\`, \`"HCR_global"\` or \`"Input"\`. Defaults to \`"HCR"\`
  under \`fmort_opt = "Catch"\` and to \`fmort_opt\` itself otherwise,
  where it is unused. Mind the indexing when mixing the two:
  \`catch_input\[r, y\]\` is the catch taken in year \`y\`, while
  \`f_ref_pt\[r, y\]\` sets F in year \`y + 1\`.

- catch_terminal_yr:

  Logical. Whether projection year 1, which replays the terminal
  assessment year, is solved against its catch target rather than fished
  at \`terminal_F\`. Default \`FALSE\`. \`TRUE\` suits the common case
  where the terminal year is not yet complete, but it overrides the F
  the assessment estimated and so changes the numbers entering year 2.
  With \`n_seas \> 1\` the terminal year takes all its seasons from
  \`terminal_NAA\`, so only the last season's F feeds year 2.

- catch_f_max:

  Upper bound on the F searched under \`fmort_opt = "Catch"\` and by
  \`F_constraint\`. Default 5. An unreachable target caps F here,
  undershoots, and warns with the regions named.

- catch_tol:

  Relative catch tolerance for the F solver. Default 1e-6.

- catch_max_iter:

  Maximum solver iterations per projection year. Default 100.

- t_spawn:

  Fraction of the spawning season elapsed before spawning.

- srr_opt:

  Named list of inputs for deterministic stock-recruit recruitment under
  \`recruitment_opt = "bh_rec"\` or \`"ricker_rec"\`, passed straight to
  [`Get_Det_Recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Det_Recruitment.md)
  and holding every argument that function needs. Formerly
  \`bh_rec_opt\`. The arrays are `R0` `[n_pop]`, `h` and
  `rec_region_prop` `[n_pop, n_regions]`, `rec_seas_prop`
  `[n_pop, n_seas]`, `SSB` `[n_pop, n_regions, n_yrs]`, `WAA`, `MatAA`
  and `natmort` `[n_pop, n_regions, n_seas, n_ages]` (natmort also
  accepted without the season dim), `Movement`
  `[n_pop, n_regions, n_regions, n_seas, n_ages]`,
  `sgl_seas_spawning_movement` `[n_pop, n_regions, n_regions, n_ages]`,
  `stray_rate` `[n_pop]`, `init_F` `[n_regions, n_seas, n_fish_fleets]`,
  `fish_sel` and `ret_sel`
  `[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]`, `dmr`
  `[n_regions, n_seas, n_fish_fleets]` and `sex_ratio_f`
  `[n_pop, n_regions]`. The scalars are `rec_dd`, `rec_lag`, `n_pop`,
  `n_regions`, `n_ages`, `n_seas`, `spawn_seas`, `seasdur`, `t_spawn`
  and `do_recruits_move`. Spawning biomass is built internally by
  appending projected SSB to `srr_opt$SSB`.

  `srr_opt$rec_lag = 1` computes each year's recruitment up front from
  the prior year's SSB, as `"inv_gauss"` and `"mean_rec"` do.
  `srr_opt$rec_lag = 0` computes it from the year's own SSB once
  `spawn_seas` is reached, and inserts recruits no earlier than that
  season, so `rec_seas_prop` must be zero before it. Reference points
  and the seasonal SBPR are unaffected either way.

- bh_rec_opt:

  Deprecated former name of `srr_opt`; supplying it warns and forwards,
  and supplying both is an error.

- n_seas:

  Integer. Number of seasons. Default 1.

- seasdur:

  Numeric vector \`\[n_seas\]\` of season durations as fractions of a
  year.

- spawn_seas:

  Integer spawning season index.

- rec_seas_prop:

  Array \`\[n_pop, n_seas\]\` of the share of annual recruitment
  entering in each season, summing to 1 for each population.

- Mrate:

  Array dimensioned like `Movement` holding the generator rather than
  the realized fractions. Only read when \`move_timing\` is 1 or 2, so
  \`NULL\` (default) is valid under \`move_timing = 0\`.

- move_timing:

  When movement happens relative to mortality within a season. \`0\`
  (default) moves then kills, \`1\` kills then moves, and \`2\` runs the
  two continuously, which also switches catch at age to the spatial
  Baranov form on season-integrated abundance. Must match the timing the
  reference points were derived under.

- expm_nsub:

  Integer controlling how the matrix exponential is evaluated under
  \`move_timing = 2\`: \`0\` uses \`Matrix::expm\`, \`n \>= 1\` uses
  \`n\` implicit backward Euler substeps. See \[mat_exp()\].

- rec_devs:

  Optional array `[n_pop, n_regions, n_proj_yrs]` of multiplicative
  deviations applied to whatever recruitment `recruitment_opt` produces,
  so a deterministic option becomes stochastic under deviations the
  caller draws. `NULL` (default) leaves recruitment as the option gives
  it. Year 1 generates no recruitment, so its slice is never read.
  Drawing outside is what lets replicates share recruitment across
  management procedures, and what lets the projection be differentiated
  with respect to the rule with the deviations kept fixed.

- naa_devs:

  Optional array
  `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]` of
  multiplicative innovations on the numbers at age, applied where a
  state-space fit applies its state: at the start of a year after ageing
  and the plus group, and at a later season's start after survival and
  movement. Unfished numbers are scaled too, as in the fit. One wherever
  there is no state. Year 1 is the terminal year replayed from
  `terminal_NAA`, so its slice is never read. A fit with
  `n_proj_yrs_devs` holds the projected innovations in `ln_NAA` past the
  fitted years, so `exp(ln_NAA[,, n_yrs + k,,,])` fills year `k + 1`.
  `NULL` (default) applies none.

- F_constraint:

  Optional constraint on the F the rule or input gives, a function
  `(F, state)` returning the annual F by region to fish at, or a list of
  them applied in the order given.
  [`catch_cap`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md)
  and
  [`catch_stability`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md)
  build the two usual ones. Applied at the front of each year from year
  2 on, since year 1 replays the terminal year, and never in a year
  solved to a catch target. `state` holds `y`, `F_rule`, `last_catch`
  and `SSB` `[n_pop, n_regions]` from the year before, `NAA`
  `[n_pop, n_regions, n_ages, n_sexes]` at the start of this year,
  `catch_at_F(F)`, the retained catch each region would take this year
  at an annual F, and `F_at_catch(target, F)`, the F that takes `target`
  where it is given and keeps `F` where it is `NA`. Catch is in the
  units of `proj_Catch`, summed over populations, seasons and fleets.
  Inside the model the F that takes the allowed catch is solved
  numerically first, then refined with two Newton steps the model
  differentiates through, so `sdreport` accounts for how that F moves
  with the estimates. A constraint of your own must make its comparisons
  through
  [`ad_value`](https://chengmatt.github.io/SPoRC/dev/reference/ad_value.md).

## Value

A named list of projected quantities. Year index 1 is the terminal
assessment year replayed, so year 2 is the first projected year and
catch advice for terminal year + 1 is read from index 2. `proj_NAA`,
`proj_NAA0`, `proj_F` and `proj_F_seas` fill their trailing
`n_proj_yrs + 1` year slot; `proj_ZAA`, `proj_ret_FAA` and
`proj_disc_FAA` leave it at 0.

- `proj_F`:

  \`\[n_regions, n_proj_yrs + 1\]\`. Annual F by region, summed over
  seasons and fleets. The trailing column holds the F the rule or input
  would apply in the year after the projection, and stays 0 under
  \`fmort_opt = "Catch"\`.

- `proj_F_rule`:

  Dimensioned as \`proj_F\`. The F the rule or input gave before any
  \`F_constraint\`, so a year where it differs from \`proj_F\` is one a
  constraint bound. 0 in years solved to a catch target.

- `proj_F_seas`:

  \`\[n_regions, n_proj_yrs + 1, n_seas\]\`. The same F by season, so
  \`rowSums(proj_F_seas\[, y, \])\` recovers \`proj_F\[, y\]\`. The only
  place the answer lives under seasonal catch targets.

- `proj_ret_FAA`:

  \`\[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes,
  n_fish_fleets\]\`. Retained fishing mortality-at-age, the component
  that lands catch.

- `proj_disc_FAA`:

  Dimensioned as \`proj_ret_FAA\`. Discard fishing mortality-at-age, set
  by \`ret_sel\` and \`dmr\`. Total F at age is the sum of the two.

- `proj_Catch`:

  \`\[n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets\]\`. Retained
  catch in biomass, built from \`proj_ret_FAA\`, and the quantity
  \`catch_input\` is matched against.

- `proj_SSB`:

  \`\[n_pop, n_regions, n_proj_yrs\]\`. Female spawning biomass,
  accumulated in \`spawn_seas\` with the \`t_spawn\` correction. Halved
  when \`n_sexes = 1\`.

- `proj_eff_SSB`:

  \`\[n_pop, n_proj_yrs\]\`. Effective spawning biomass at each
  population's natal region, with cross-population terms scaled by
  \`stray_rate\`. Equal to SSB summed across regions when \`n_pop = 1\`.

- `proj_Total_Biom`:

  \`\[n_pop, n_regions, n_proj_yrs\]\`. Total biomass over all ages and
  sexes, at the same point in the season as \`proj_SSB\` and on the
  estimation model's definition, so the series continues without a
  discontinuity at the terminal year.

- `proj_Dynamic_SSB0`:

  \`\[n_pop, n_regions, n_proj_yrs\]\`. Spawning biomass under the same
  realized recruitment but no fishing, for dynamic depletion.

- `proj_NAA`:

  \`\[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes\]\`.
  Fished numbers at age at the start of each season, before that
  season's mortality and ageing. Movement has already been applied under
  \`move_timing = 0\` and not under 1 or 2. The trailing slot holds the
  numbers passed into the year after the projection.

- `proj_NAA0`:

  Dimensioned as \`proj_NAA\`, decremented by natural mortality alone.

- `proj_ZAA`:

  Dimensioned as \`proj_NAA\`. Total mortality-at-age for the season:
  natural mortality scaled by season duration plus retained and discard
  F summed over fleets.

- `proj_catch_resid`:

  Shaped like \`catch_input\`. Relative miss on each target,
  \`(realized - target) / target\`, and \`NA\` for years with no target.
  Should be at or below \`catch_tol\` wherever the solve converged, and
  is worth checking directly rather than relying on warnings.

## Details

A projection year generates and allocates recruitment, builds F at age
from the annual F, the terminal year's seasonal ratios and selectivity,
moves fish each season, applies within-season mortality and ages the
survivors at the end of the final season, computes spawning biomass in
`spawn_seas` with a mid-season correction (after spawning movement under
natal homing with one season), takes catch by the Baranov equation, and
sets next year's F from the control rule or input.

Under `srr_opt$rec_lag == 0` the recruitment and spawning steps are
reordered within `spawn_seas`: movement runs first, spawning biomass is
computed from the survivors alone, that SSB generates this year's
recruitment, and the recruits are inserted immediately before mortality
and ageing. Year 1 advances the terminal state with no recruitment
event.

Under `fmort_opt = "Catch"` the F step moves to the front of the
following year, since the F that lands a target depends on that year's
own numbers at age rather than the previous year's spawning biomass. The
year is run at trial F values until realized catch matches the target,
then run once more at the accepted F and committed; no demographic input
is modified. Regions are solved jointly, because movement (and, under
`move_timing = 2`, the season-integrated abundance) makes each region's
catch depend on every other region's F. Seasonal targets are swept
forward one season at a time, which is exact because a season's catch
depends only on the F in that season and earlier ones. The catch solved
against is retained catch, so a fleet that discards exerts more total F
than the target implies.

Under `F_constraint` the F the rule or input gave is checked at the
front of the year in the same way, from year 2 on. Each constraint reads
the retained catch that F would take this year, through the same trial
runs, and may lower or raise F by solving for the catch it allows.
[`catch_cap`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md)
and
[`catch_stability`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md)
build the two usual ones.

Spawning biomass is multiplied by 0.5 when `n_sexes = 1`, and movement
is skipped when `n_regions = 1`.

## Differentiating through the projection

The projection uses RTMB's replacement operators, so it can run inside
[`MakeTape`](https://rdrr.io/pkg/RTMB/man/Tape.html) or `MakeADFun` and
be optimized with an exact gradient, giving an F schedule solved against
an objective rather than scanned over a grid.

Catch targets (`fmort_opt = "Catch"`) and `F_constraint` run inside the
model too. The F that takes a given or an allowed catch is solved on the
plain numbers, then refined with two Newton steps the model
differentiates through, so the F has the derivative of the solve. RTMB
records the model once, when it is built, and keeps whichever side of
each comparison was true then, so a model that must compare at the
estimates is built again there before `sdreport`, as
[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
does.

One option is refused inside the model, since it has no derivative and
would otherwise return a wrong gradient rather than an error:
`recruitment_opt = "inv_gauss"` draws at random. Use `"mean_rec"`,
`"bh_rec"` or `"ricker_rec"` there. A control rule that compares
spawning biomass directly stops with an error inside RTMB. One that
reads it through
[`ad_value`](https://chengmatt.github.io/SPoRC/dev/reference/ad_value.md),
as
[`hcr_threshold`](https://chengmatt.github.io/SPoRC/dev/reference/hcr_threshold.md)
does, runs, with the side of the comparison fixed where the model was
built.

## See also

Other Reference Points and Projections:
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`Setup_Mod_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md),
[`catch_cap()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md),
[`catch_stability()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md),
[`hcr_threshold()`](https://chengmatt.github.io/SPoRC/dev/reference/hcr_threshold.md)
