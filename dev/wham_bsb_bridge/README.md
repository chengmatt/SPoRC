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
| Fleets | commercial and recreational, in each region, kept as four fleets so each takes its own composition likelihood |
| Surveys | recreational CPA and VAST, in each region, also four |
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
| `06_figures.R` | writes the bridge figures |
| `07_refit.R` | builds the estimation model, initialized the way wham initializes |
| `08_fit.R` | optimizes it |
| `09_compare_fits.R` | the usual SPoRC plots, refit against wham |
| `10_q_fixed_check.R` | refits with catchability held at wham's, the diagnostic for the northern gap |

The scripts load SPoRC from source, since the seasonal reporting option this bridge
uses is newer than any installed build.

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

## Two likelihoods that had to be added

WHAM fits one annual catch figure and one annual age composition per fleet, summed
over all 11 seasons, while the survey indices are snapshots inside one season.
SPoRC used to fit every observation in the season it sat in, so the fishery
likelihood could not be written against the same observations at all. Six of the
eight compositions are also on WHAM's logistic normal with the empty bins dropped,
which SPoRC did not have. Both are now options.

## Reporting once a year

Every data source now takes a `_seas_Type`: `"spltSeas"` fits the season, `"aggSeas"` sums the
prediction over the year and fits one observation. This bridge sets the catch and
the fishery compositions to `"aggSeas"` with the observation in season 1, and
leaves the survey indices and their compositions seasonal.

## Figures

`figures/` holds six, all with every parameter held at the WHAM values:

| File | What it shows |
|---|---|
| `01_ssb_recruitment.png` | spawning biomass and recruitment by stock, both models |
| `02_predicted_catch.png` | predicted catch by fleet |
| `03_survey_indices.png` | predicted survey indices with the observations behind them |
| `04_index_difference.png` | the index difference as a percent, zero outside season 4 |
| `05_relative_error.png` | how closely each quantity bridges, on a log scale |
| `06_likelihood_components.png` | negative log likelihood by component |

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

### Likelihood components

Negative log likelihood by fleet, in `output/05_likelihood_components.csv`:

| Component | Fleets that match | Fleets that do not |
|---|---|---|
| Aggregate catch | all four, to 2e-12 | none |
| Survey index | 1 and 3, to 1e-12 | 2 and 4, by the movement ordering below |
| Fishery age compositions | all four, to 5e-10 | none |
| Survey age compositions | 1 and 3, to 3e-09 | 2 and 4, by the movement ordering below |

Two adjustments got the compositions there. Fleet 1's fishery compositions and
survey 2's are on a Dirichlet multinomial, where WHAM uses the saturating form,
concentration \(\alpha_a = p_a \exp(\theta)\), and SPoRC the linear form,
\(\alpha_a = p_a \exp(\theta) N\); the two agree once the input sample size is
divided out of the dispersion, which is what `04_seed_at_wham.R` does. The other six
are on WHAM's logistic normal with the empty bins dropped, which SPoRC now offers as
`"iid-Logistic-Normal-miss0"` and `"1d-Logistic-Normal-miss0"`; that form takes
WHAM's parameters directly, with no adjustment.

Every remaining difference is on a survey sitting in season 4, which is the
movement ordering below and nothing else.

## Dropping the empty bins

WHAM's `logistic-normal-miss0` and `logistic-normal-ar1-miss0` drop the empty bins,
renormalize the expected proportions over the bins that remain, divide the standard
deviation by the square root of the input sample size, and subtract the change of
variables from the log ratio so the result is a density on the composition. SPoRC's
existing logistic normal forms add a constant to the empty bins and do none of the
other three, so they are different functions of the same data rather than the same
function parameterized differently. The two new options reproduce WHAM's form and
take its parameters unchanged.

## Refitting in SPoRC

`07` builds the estimation model. It mirrors wham where SPoRC can: one annual fishing
mortality per fleet, the eleven seasons of a year sharing a deviation; four survey
catchabilities; mean recruitment and a recruitment deviation scale per stock; a
numbers at age deviation scale estimated for each stock in its own region and fixed
at 0.05 for north fish sitting in the south; the 2dar1 age and year correlations, on
SPoRC's scale two transform where wham uses scale one; the composition dispersion and
correlation per fleet; and an equilibrium initial age structure from one recruitment
level and one initial fishing mortality per stock, which is wham's `N1_model = 1`.
Movement, selectivity, natural mortality and the observation errors stay at wham's.

162 fixed effects and 962 random effects. It converges to a maximum gradient of
1.1e-10 with a positive definite Hessian.

Recruitment deviations follow the same first order autoregression wham gives them, and
the first year is estimated without a penalty, since it is the first year's age one
abundance and belongs to the initial condition. Two things wham does that SPoRC does
not: selectivity is fixed here rather than estimated with its own deviations, and the
environmental covariate is absent, which costs nothing in `fit_0`.

### What the refit gives

| | North | South |
|---|---|---|
| Spawning biomass, mean absolute difference | 8.8% | 0.7% |
| Spawning biomass, terminal | +8.8% | -0.6% |
| Recruitment, mean absolute difference | 8.9% | 1.4% |

The southern stock lands on wham. The northern stock sits about 8% above it, and the
gap is a catchability difference and nothing else: the two northern catchabilities are
estimated 8.6% and 6.7% below wham's while the two southern ones match to 0.4% and
0.01%. Catchability down and biomass up in the same proportion leaves the index fits
identical, which is why the index and catch figures show one curve.

Holding catchability at wham's values and refitting everything else closes it. The
northern mean absolute difference falls from 8.8% to 1.7% and the terminal year from
+8.8% to -0.4%, while the objective moves by 0.3 units. The northern scale is close to
flat in this likelihood, so the two models settle at different points along it without
either fitting the data better.

Two explanations that turned out not to be the cause. Adding the recruitment
correlation moved the northern gap from 8.5% to 8.8%, so it is not that, although it
did recover wham's correlation and brought the northern recruitment scale from -0.020
to -0.102 against wham's -0.080. And the movement ordering below is worth 0.2 to 0.5%
on the two affected indices at the refit's own numbers, far too little to carry 8%.

Parameters, SPoRC refit against wham:

| Parameter | SPoRC refit | WHAM |
|---|---|---|
| log deviation scale, north / south | -0.199 / -0.168 | -0.215 / -0.168 |
| log recruitment scale, north / south | -0.102 / -0.597 | -0.080 / -0.599 |
| recruitment correlation, north / south | 0.376 / -0.044 | 0.361 / -0.016 |
| age correlation, north / south | 0.098 / 0.252 | 0.137 / 0.240 |
| year correlation, north / south | 0.235 / 0.128 | 0.278 / 0.121 |
| log catchability, four surveys | -9.22, -4.17, -8.84, -4.19 | -9.13, -4.10, -8.85, -4.19 |

The two joint likelihoods are not comparable, since wham's holds the environmental
covariate and SPoRC's states the recruitment penalty differently.


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

The environmental covariate is not in SPoRC at all. In `fit_0` the temperature has no link to
recruitment or natural mortality, so it contributes an additively separable likelihood and no
population dynamics, and leaving it out changes nothing here.
