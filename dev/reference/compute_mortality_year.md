# Fishing and total mortality for one year

Also derives the selection-weighted weight at age where a fleet asks for
it. Both are done a year at a time because under cohort growth the key a
length-based selectivity acts through is only known once the population
loop reaches that year; every other model runs them for all years before
the loop starts. Under cohort growth this runs inside the population
loop, so the year's state is taken as an argument and handed back rather
than assigned into this frame.

## Usage

``` r
compute_mortality_year(
  y,
  st,
  growth_model,
  derive_waa,
  fish_selex_type,
  ret_selex_type,
  srv_selex_type,
  fish_waa_selected,
  srv_waa_selected,
  fish_sel_l,
  ret_sel_l,
  srv_sel_l,
  wt_len_pars,
  growth_len_mid_vals,
  UseCatch,
  UseCatch_pop,
  missing_catch,
  ln_F_mean,
  ln_F_devs,
  logit_dmr_mean,
  logit_dmr_devs,
  SizeAgeTrans,
  natmort,
  seasdur,
  n_pop,
  n_regions,
  n_seas,
  n_ages,
  n_sexes,
  n_fish_fleets
)
```

## Arguments

- y:

  Year index.

- st:

  Named list carrying `Fmort`, `dmr`, `fish_sel`, `ret_sel`, `ret_FAA`,
  `disc_FAA`, `tot_FAA`, `ZAA`, `WAA_fish`, `WAA_srv`,
  `SizeAgeTrans_fish` and `SizeAgeTrans_srv`, returned with year `y`
  updated. The last two may be `NULL`, in which case the shared
  `SizeAgeTrans` key is used instead.

- growth_model, derive_waa:

  Growth module switches.

- fish_selex_type, ret_selex_type, srv_selex_type:

  Integer (0/1); `1` for length-based selectivity.

- fish_waa_selected, srv_waa_selected:

  Integer vectors `[fleet]` (0/1) for which fleets get a
  selection-weighted weight at age.

- fish_sel_l, ret_sel_l, srv_sel_l:

  Arrays `[region, year, len, sex, fleet]` of selectivity at length.

- wt_len_pars:

  Array `[pop, region, sex, 2]` of weight-length parameters.

- growth_len_mid_vals:

  Length bin midpoints the weight-length relationship is read at.

- UseCatch, UseCatch_pop:

  Arrays flagging which cells fit an aggregate/pop-specific catch
  observation.

- missing_catch:

  Logical array, `TRUE` where the aggregate catch observation is missing
  (not a true recorded zero).

- ln_F_mean, ln_F_devs:

  Log fishing mortality mean and deviations.

- logit_dmr_mean, logit_dmr_devs:

  Logit discard mortality rate mean and deviations.

- SizeAgeTrans:

  Shared size-age key, used when no growth-derived per-fleet key is
  supplied in `st`.

- natmort, seasdur:

  Natural mortality at age and season duration.

- n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets:

  Dimensions.

## Value

`st` with year `y` updated.
