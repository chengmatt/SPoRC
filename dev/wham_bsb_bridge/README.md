# Multi-WHAM black sea bass bridge

Bridges `fit_0` from Miller, Curti and Hansell (multi-region, multi-stock WHAM applied to
northern and southern black sea bass) into SPoRC. `fit_0` is the base run: no temperature
effect on recruitment or natural mortality, no natural mortality random effects.

## What the model is

| | |
|---|---|
| Years | 1989-2021 (33) |
| Ages | 1-8, plus group at 8 |
| Populations | BSB North (natal region North), BSB South (natal region South) |
| Regions | North, South |
| Seasons | 11, one month each apart from a two month June-July step |
| Sexes | 1 |
| Fleets | commercial and recreational, in each region |
| Surveys | recreational CPA and VAST, in each region |
| Movement | only the north stock moves; north to south in seasons 7-11, south to north in seasons 1-4, everything forced back north at the end of season 5 |
| Spawning | halfway through season 6, so 0.5 of a year |
| Numbers at age | state-space, 2dar1 across ages and years, recruitment decoupled |

## Scripts

Run in order. Each writes to `output/`.

| Script | What it does |
|---|---|
| `01_wham_fit0.R` | rebuilds the WHAM input and pastes the published parameter list back in |
| `02_wham_to_sporc_data.R` | reshapes the WHAM data and report into SPoRC arrays |
| `03_build_sporc.R` | builds the SPoRC input list |
| `04_seed_at_wham.R` | holds every SPoRC parameter at the WHAM values and reports |
| `05_compare.R` | writes the comparison tables |

## Reproducing the WHAM fit

The paper repository ships parameter lists but no fitted objects, so `01` rebuilds the input
and evaluates at `parLists_no_M_re.RDS[[1]]` with `do.fit = FALSE`. That lands on
nll = -975.0993 with a maximum gradient of 5.7e-11, so it is the published maximum
likelihood estimate and not a nearby point.

WHAM 2.1.0.9009 and 9010 both have a bug that this model trips: `selectivity$n_selblocks`
overrides the ASAP block count, but `set_selectivity` leaves `selpars_ini` at the ASAP row
count, so comparing it against the bounds gives non-conformable arrays. The BSB fits collapse
12 ASAP blocks to 8, so `01` trims the stale rows and reinstalls the function.

## How the two models line up

WHAM runs four fleets and four indices as region by gear. SPoRC runs gear across regions, so
two fleets and two surveys cover the same four of each. Selectivity, weight at age and
catchability are all region-specific in SPoRC, so nothing is lost.

WHAM's `mig_type = 0` is survival over the season and then instantaneous movement, which is
SPoRC's `move_timing = 1`. The transition matrices are taken straight out of the WHAM report
through `use_fixed_movement`, so the `must_move` and `can_move` structure comes along without
having to be respecified.

`Fmort` in SPoRC is already integrated over the season, while natural mortality is an annual
rate scaled by season duration. WHAM's F is an annual rate throughout, so
`Fmort = F_wham * seasdur`. Survey timing and spawning timing are fractions of their season in
SPoRC and fractions of a year in WHAM, so both are divided by the season duration.

Year one comes from `init_age_strc = "free"`, which makes the initial deviations the log
numbers at age themselves, and the state-space cells for years two onward are set to WHAM's
realized numbers. Recruitment deviations are calibrated against a first pass so age one lands
on WHAM's numbers exactly. Everything is then held, so this run tests the seasonal spatial
transition algebra rather than any fitting.

## Result

Maximum relative error against the WHAM report, in `output/05_bridge_comparison.csv`:

| Quantity | Cells | Max relative error |
|---|---|---|
| Numbers at age on January 1 | 759 | 9.0e-16 |
| One step ahead predicted numbers at age, ages 2+ | 672 | 2.1e-15 |
| Spawning biomass | 66 | 1.5e-13 |
| Predicted catch, four fleets | 132 | 1.1e-12 |
| Predicted catch at age, four fleets | 1056 | 6.9e-11 |
| Predicted index, recreational CPA north and south | 66 | 2.5e-13 |
| Predicted index, VAST north and south | 66 | 1.1e-02 |
| Predicted index, VAST, WHAM movement order | 66 | 4.9e-13 |

## The one convention gap

The two VAST indices sit in season 4, which is a movement season, and they are the only
quantities that do not bridge. WHAM survives the fish for the index timing and then applies
the whole season's movement matrix before reading the index, so the index is measured at the
post-movement location. SPoRC moves at the end of the season, so at a timing of zero the index
is measured before the move. Reordering SPoRC's numbers the WHAM way brings both indices back
to 4.9e-13, which is in `05_compare.R`.

The recreational CPA indices are in season 6, where no movement is allowed, so both orderings
agree and they bridge without any adjustment.

## What is not bridged

WHAM fits one annual catch observation and one annual age composition per fleet, summed over
all 11 seasons. SPoRC fits catch and compositions per season and has no season-aggregated
option, so the fishery likelihood cannot be written down against the same observations. The
survey indices and survey age compositions each sit in one season and do map across.

The environmental covariate is not in SPoRC at all. In `fit_0` the temperature has no link to
recruitment or natural mortality, so it contributes an additively separable likelihood and no
population dynamics, and leaving it out changes nothing here.
