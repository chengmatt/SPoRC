# GOA rex sole bridge

SPoRC rebuilt from the 2025 Gulf of Alaska rex sole assessment (Stock Synthesis
Model 25.1, `runs/2025_models/run_14_asfor7_ageing_error` in Carey McGilliard's
`goa_rex` repository) and evaluated at the assessment's own maximum likelihood
estimate. Two areas with growth estimated separately in each, survey conditional
age-at-length, an ageing error matrix.

## Files

- `R/build_rex_data.R` parses the SS3 input files with r4ss into SPoRC arrays
  (`build_rex_data()`), and attaches the report quantities and parameter estimates
  from `Report.sso` (`add_rex_ss3_report()`). Writes `output/rex_bridge_data.rds`.
  The packaged copy is `data/mlt_rg_goa_rex_data.rda` (raw parsed files and the
  bulky per-row tables dropped).
- `R/run_rex_bridge.R` builds the SPoRC model through
  `tests/testthat/helper-bridge_goa_rex.R`, seeds the estimate, and prints every
  comparison: growth tables, age-length key, selectivity, numbers at age, spawning
  biomass, recruitment, total biomass, catch, indices, likelihood components and the
  gradient at the estimate.
- `tests/testthat/test-regression_goa_rex_bridge.R` is the regression gate.

## Running SS3

The run folders in the assessment repository hold inputs and a Windows
executable only. `Report.sso` was produced by compiling SS3 v3.30.18 from
`nmfs-ost/ss3-source-code` with ADMB 12.3 (`./Make_SS_330_new.sh -a <admb> -b
build -o`) and running `ss_opt -nohess` in a copy of the run folder. The current
SS3 main needs ADMB 13; built with 12.3 its raw data source reads fail on the data
file (a garbage size-frequency count, then an out-of-memory kill).

## What matched, and what does not

At the assessment's estimate: numbers at age, spawning biomass, recruitment, total
biomass, catch and both indices to 0.005% in both areas (the assessment's
six-digit parameter rounding); mean length, SD and weight at age to 5e-4%; the
age-length key to 2e-6; every selectivity curve to 3e-6. Composition likelihoods
agree up to the factor (1 + n_bins * addtocomp) from SS3 renormalizing after
adding its constant; the index likelihood up to 0.5 log(2 pi) per observation;
the catchability prior up to log(sd) + 0.5 log(2 pi).

Not reproduced, by convention: the recruitment and early-deviation penalties
(SS3 applies each early year's own bias adjustment, SPoRC centers the initial
deviations on the first ramp value), (cross-fleet catchability mirroring is expressed through the map, which ties the
two surveys' q cells to one parameter).
Those leave the gradient at the estimate non-zero on R0, the growth parameters and
the apportionment, without moving any of the quantities above.

## Refit conventions

At the SS3 estimate every likelihood component reproduces (see the regression
test). On refit the two part company only through `do_recdev = 1`: SS3's main
deviations are an ADMB dev_vector constrained to sum to zero, SPoRC's are free
under the penalty, so R0 and the deviation mean trade off in SPoRC (R0 about
4% lower, mean deviation +0.03). A finite-difference check of SS3 against
SPoRC (`scratchpad/fd_run.sh`, `fd_sporc.R`, `fd_collect.R`; SS3 re-evaluated
at perturbed par files with `-maxfn 0`) found every other parameter's
likelihood partials identical to print precision.

Checked directly on 2026-08-23: a throwaway build imposing the sum-to-zero
constraint (last main deviation determined, still penalized) refits to SS3
within 0.02% SSB / 0.04% recruitment, ln R0 11.53128 vs 11.53150, deviations
within 2e-4, objective 0.005 below the SS3 estimate. The constraint is not kept
in the package.
With the catch CV tightened to 1e-4 as well (hybrid F matches catch exactly), the
same build agrees to 0.0004% SSB / 0.001% recruitment; renormalizing the
composition constant (weights 1/(1 + n c)) moves only 0.001%.
