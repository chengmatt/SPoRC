# EBS Pacific cod bridge

SPoRC rebuilt from the 2024 eastern Bering Sea Pacific cod assessment (Stock
Synthesis Model 24.1, Barbeaux et al. 2024) and evaluated at the assessment's own
maximum likelihood estimate. One area, one sex, ages 0-20, 1977-2024, with
time-varying growth carried cohort by cohort, length-based double normal
selectivity on 121 population length bins fit to 24 five-centimetre data bins,
and annual deviations on the survey's ascending width.

This is the bridge that motivated the time-varying growth module: the assessment
varies the length at age 1.5 and the Richards K parameter annually from 2000,
each on the bounded logit scale, and propagates size at age through cohorts
rather than reading it off each year's curve.

## Files

- `R/build_pcod_data.R` parses the SS3 input files with r4ss into SPoRC arrays
  (`build_pcod_data()`), and attaches the report quantities, the year-by-year
  weight and fecundity at age from `wtatage.ss_new`, and the parameter estimates
  from `ss.par` (`add_pcod_ss3_report()`). Writes `output/pcod_bridge_data.rds`;
  the packaged copy is `data/sgl_rg_ebs_pcod_data.rda`.
- `tests/testthat/helper-bridge_ebs_pcod.R` builds the SPoRC input list
  (`build_ebs_pcod_input()`) and seeds every parameter at the assessment's
  estimate (`seed_ebs_pcod_mle()`).
- `tests/testthat/test-regression_ebs_pcod_bridge.R` is the regression gate,
  57 assertions against the assessment's own reported quantities.
- `vignettes/ae_ebs_pacific_cod_case_study.Rmd` is the walkthrough: one section
  per `Setup_Mod_*` call in the order the helper makes them, with the reason for
  each argument that is not a `SPoRC` default. Its code builds a model identical
  to the helper's, checked slot by slot.
- `dev/make_sporc_obj_figs/make_ebs_pcod_bridge_figs.R` builds the model through
  the same helper, refits it, and writes the case study figures to
  `vignettes/figures/ae_ebs_pcod_*.png`: the time series at the assessment's
  estimate and refitted, selectivity at length for both fishery blocks and four
  survey years, mean length at age across the time-varying period, the two
  varying growth parameters through the series, and expected length
  compositions for both fleets.

## Running SS3

The assessment repository ships inputs only. `Report.sso` and `wtatage.ss_new`
were produced by compiling SS3 v3.30.18 from `nmfs-ost/ss3-source-code` with
ADMB 12.3 (`./Make_SS_330_new.sh -a <admb> -b build -o`) and running
`ss_opt -maxfn 0 -nohess` in a copy of the Model 24.1 folder with
`init_values_src` set to 1, so the run evaluates at `ss.par` without estimating.
That reproduces the assessment's own reported objective exactly: 174 parameters,
total 246.832, with every component matching the assessment document's Table
2.14.

## What matched

At the assessment's estimate, every comparison is at the report's six-digit
print precision:

| Quantity | Agreement |
|---|---|
| Weight at age (population, fishery-selected, mid-year) | 5e-4 % |
| Growth parameters by year | 1e-3 % |
| Selectivity at length, both fleets, every year | 5e-7 |
| Selectivity folded to age, both fleets, every year | 9e-7 |
| Numbers at age | 0.001 % |
| Spawning biomass, recruitment, total biomass, catch | 5e-4 % |
| Survey index | 4e-4 % |
| Expected length and age compositions | 1e-5 |

and the likelihood components, each restated in the assessment's convention:

| Component | SPoRC | SS3 |
|---|---|---|
| Survey index | -59.065728 | -59.0658 |
| Fishery lengths | 119.364878 | 119.365 |
| Survey lengths | 103.268624 | 103.269 |
| Survey ages | 57.440030 | 57.4390 |
| Recruitment (main and early) | 5.0657970 | 5.06581 |
| Initial equilibrium recruitment offset | 1.9070368 | 1.90704 |
| Growth deviations | 10.415703 | 10.4157 |
| Survey selectivity deviations | 8.4322459 | 8.43225 |

The catch is fit essentially exactly on both sides (SS3 6.9e-11, SPoRC 5.4e-7
after removing the normal constants), since the assessment conditions fishing
mortality on the catch.

## Parameter counts

The assessment reports 174 active parameters and SPoRC estimates 220. Every
block corresponds one for one except fishing mortality, initial fishing
mortality and the survey's extra standard deviation.

| Block | SS3 | SPoRC |
|---|---|---|
| Virgin recruitment | 1 | `ln_global_R0` 1 |
| Initial equilibrium offset | `SR_regime` 1 | `ln_rinit` 1 |
| Survey catchability | 1 | `ln_srv_q` 1 |
| Survey extra SD | 1 | folded into the index standard error, 0 |
| Initial fishing mortality | 1 | held at the assessment's value, 0 |
| Growth curve (`L1`, `L2`, `K`, Richards) | 4 | `ln_growth_pars` 4 |
| Growth deviations (`L1` and `K`, 2000-2024) | 50 | `ln_growth_devs` 50 |
| Survey ascending-width deviations | 43 | `ln_srvsel_devs` 43 |
| Selectivity (4 fishery, 2 survey) | 6 | 4 + 2 |
| Main recruitment deviations | 46 | `ln_RecDevs` 46 |
| Early (initial age) deviations | 20 | `ln_InitDevs` 20 |
| Fishing mortality | 0, solved from catch | `ln_F_devs` 48 |
| **Total** | **174** | **220** |

`174 + 48 - 1 - 1 = 220`. The 48 are the whole difference in kind: the
assessment conditions fishing mortality on the catch with its hybrid solver, so
F costs it no parameters, while SPoRC estimates a deviation per year against a
soft catch likelihood. The other two are held rather than dropped.

## Refit

Optimizing from the assessment's estimate converges to a maximum gradient of
1.3e-9 with a positive definite Hessian (minimum eigenvalue 0.033), and drops
the objective 0.73. Spawning biomass lands a median +0.55
percent from the assessment (max 3.2) and recruitment +0.19 percent. Growth
comes back within 0.1 percent on every parameter, which is the check that
matters for this bridge: `L1` 13.836 against 13.851, `L2` 112.23 against 112.26,
`K` 0.11431 against 0.11453, the Richards coefficient 1.4881 against 1.4856.

**The plus group's initial deviation needs `equil_init_age_strc = "stoch_all"`.**
The assessment estimates and penalizes all twenty early deviations, the one at
the accumulator age included. SPoRC's `"stoch_plus_grp"` instead holds that
deviation at zero and leaves it out of the initial-age penalty, so the plus group
follows the equilibrium calculation. Both are coherent settings, and SPoRC keeps
the map, `data$map_ln_InitDevs` and `init_devs_pen_use` in step with whichever is
chosen: a deviation is penalized exactly when it is estimated.

The trap is overriding `map$ln_InitDevs` after setup to estimate the deviation
while the two data-side copies built from the setup map are left behind. That
estimates it without penalizing it, and it is then a parameter with neither a
penalty nor any data behind it: the fish it scales were age 20 in 1977 and are
long gone by the end of the series. Its curvature measures 5e-7 against 2.1
(= 1/sigmaR^2) for every other initial deviation, it alone takes the Hessian's
smallest eigenvalue to 2e-9 and its condition number to 4e15, and the optimizer
stops with a stale gradient of 0.029 while `solve()` calls the Hessian singular.
Choosing the option at setup rather than patching the map afterwards leaves all
twenty estimated and penalized, the smallest initial-deviation curvature at 2.16,
and the refit converging to 1.3e-9.

The drift is the recruitment convention, not a specification error. The
assessment declares `do_recdev = 1`, an ADMB `dev_vector` constrained to sum to
zero, so its R0 is the geometric mean recruitment by construction (its main
deviations sum to -5e-14). SPoRC's deviations are free under the penalty, so R0
and their mean slide together: on refit the mean deviation is -0.110 and ln R0
falls 0.106, from 13.341 to 13.235, leaving their product essentially unchanged.
That is the same convention difference recorded on the GOA rex sole and BSAI
northern rockfish bridges, and it is deliberate -- dev_vectors are a legacy
identifiability device rather than statistical practice the package carries.

The gradient at the assessment's estimate is large on the growth parameters
(11.6), the recruitment deviations (10.4) and the fishing mortality deviations
(5.4) although every reported quantity matches. That is the hybrid F: the
assessment conditions fishing mortality on the catch rather than estimating it,
so SPoRC's partials omit the implicit dF/dtheta and its soft catch likelihood
stands in for the exact solve.

## SS3 conventions that fail silently

Each of these was found by a component that would not reproduce, and each is
worth carrying to the next Stock Synthesis bridge.

- **A fishing fleet's compositions are assigned to MID SEASON whatever month the
  data file records.** `get_data_timing()` in `SS_global.tpl` sends any fleet
  with `surveytiming < 0` and no explicit timing down an "assign to midseason"
  branch that overwrites both the sub-season index and the timing fraction with
  the mid-season ones. The fishery length compositions here carry month 1, and
  reading them at the start of the year (as the month says) leaves the fishery
  length likelihood at 198.6 against the assessment's 119.4. On the mid-season
  key it is 119.3649. The numbers behind them are still the season-long Baranov
  average, not numbers decayed to mid season.
- **Length selectivity is applied AT LENGTH, not folded to age first.** SS3 builds
  `exp_AL(a, l) = N_a Z_a^{avg} P(l|a) s(l)` and sums over ages, which keeps the
  covariance of length and selection within an age. Expanding the catch at age
  through the key instead (SPoRC's default) is a different quantity. This is
  `FishLenComps_sel = "length"` / `SrvLenComps_sel = "length"`.
- **The catch in biomass uses the SELECTION-weighted weight at age**,
  `fish_body_wt = (ALK · s(l) · w(l)) / (ALK · s(l))`, not the population mean
  weight. This is `fish_waa_selected = 1`.
- **The CV at age is frozen at the START YEAR and never updated.** `Make_AgeLength_Key`
  computes `CV_G` only under `y == styr`; from then on the spread at age is that
  frozen CV times the current year's mean length. Interpolating the CV on the
  current year's own curve instead moves weight at age by up to 1.4 % in the
  time-varying years and the catch by 0.44 %, which is what localized it. Both
  the reference sizes and the `L1`/`L2` endpoints of the interpolation must come
  from the start year.
- **Time-varying growth is propagated cohort by cohort, not read off each year's
  curve.** Every cohort grows from the size it reached by the increment the
  current year's parameters imply; ages still in the linear phase keep the length
  at `A1` their birth year's parameters gave them; the first integer age past
  `A1` sits on the current year's curve; and the plus group's size next year is
  the numbers-weighted blend of the cohort entering it with the fish already
  there, `((N_{n-1} + 0.01) g(L_{n-1}) + (N_n + 0.01) g(L_+)) / (N_{n-1} + N_n + 0.02)`.
  That blend is why the module has to run inside the year loop.
- **A bounded parameter's deviations are logit deviations inside its bounds**,
  `P_y = lo + (hi - lo) logit^{-1}(logit((P - lo)/(hi - lo)) + delta_y)`, not
  multipliers. The deviation itself is unit normal and the control file's
  `dev_se` scales it.
- **`exploitation$annual_F` is not the fishing mortality.** It is an
  exploitation-rate summary; the rate that multiplies selectivity is the
  fleet-named column of the same table (equal to `F_std` under
  `F_report_units = 3`). Seeding the wrong one is a 0.5 % error in F that
  compounds to 6 % in the oldest ages by the end of the series.
- **The double normal's ascending limb is anchored at the first DATA bin**, not
  the first population bin, and bins below it take `(b / b_start)^2` times the
  selectivity there. This is `fish_sel_dbnrml_startbin`.
- **`Biology_at_age_in_endyr`'s `Len:_f` columns are the selection-weighted mean
  length of fleet `f`**, not the population mean; and r4ss's `growthseries` keeps
  only sub-season 1, so comparing a mid-season quantity against it looks like a
  7 cm error that is not there. `wtatage.ss_new` is the reliable year-by-year
  source: fleet 0 is the start-of-season population weight, -1 the mid-season
  weight, -2 the fecundity, and a positive fleet its selection-weighted weight.

## Convention differences, not reproduced

- **Single-sex spawning biomass.** SPoRC's is the female share at an even sex
  ratio; a one-sex SS3 model's spawning output counts every mature fish. The two
  differ by exactly two, and since steepness is one here spawning biomass enters
  no likelihood, so nothing else moves. The bridge compares `2 * SSB`.
- **Normal constants in the deviation penalties.** SPoRC carries
  `log(sigma) + 0.5 log(2 pi)` per deviation; SS3 carries only the bias-ramp
  share of `log(sigma)`. With every deviation standard deviation fixed, as here,
  that is a constant.
- **The bias ramp on `log(sigmaR)`.** SPoRC applies half the ramp to the
  `log(sigmaR)` term of the recruitment penalty where SS3 applies the whole ramp
  (`recr_like = sd_offset_rec * log(sigmaR) + sum(dev^2)/(2 sigmaR^2)` with
  `sd_offset_rec` the summed bias adjuster). A constant while sigmaR is fixed;
  the two would disagree materially if it were estimated.
- **The survey's extra standard deviation is added, not estimated.** SS3 fits it
  as a free parameter, so how much weight the index carries is decided with
  everything else; here it is added to `ObsSrvIdx_SE` and held. Identical at the
  assessment's estimate, since the total standard error is the same and the index
  likelihood reproduces to six digits, but on refitting SS3 can rescale the
  index's influence and the bridge cannot. One of the two parameters held rather
  than estimated in the count above.
- **Equilibrium catch (1.7e-4) and soft bounds (4.2e-3).** SPoRC has neither;
  initial fishing mortality is held at the assessment's estimate rather than
  fit to the equilibrium catch.
