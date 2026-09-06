> **Superseded in part.** The `srv_q_key` / `fish_q_key` catchability-at-age
> parameters described below were removed. A free catchability per age and a
> selectivity estimated at age are one quantity written two ways, so the age
> shape now lives in selectivity, through the `"nonparfree"` form: non-parametric
> on the log scale with no standardization. Grouping that a key expressed is done
> by `*_sel_nonpar_est_bins`. Everything else here still holds.

# Age-disaggregated observations, age-specific catchability, and effort-driven F

Design note. Nothing here is implemented yet.

## What this is for

`SPoRC` fits catch as an aggregate plus a composition, and an index as an
aggregate plus a composition. The SMS family, and SAM, do neither: they fit
catch at age and index at age directly, every age its own observation, with its
own catchability and its own variance. That is the native form for the four
North Sea sandeel stocks, for sprat, and for the SAM lineage generally.

The two are not interchangeable. The exact factorization of an at-age
observation into a total and a composition holds for Poisson and multinomial,
where the total is a genuine count sum. Both `smsR` and SAM are lognormal, and
the sum of lognormals is not lognormal, so no choice of composition likelihood
recovers the at-age statement. Measured on the North Sea sandeel bridge, the
direct formulation reproduces the reference to 6% on spawning biomass while
every aggregate variant sits between 65% and 140%, whether the compositions are
multinomial, Dirichlet-multinomial or logistic-normal. That gap is structural.

Today the only way to express it is to stand up one fleet per age with
selectivity pinned, which is what
`vignette("ah_north_sea_sandeel_case_study")` does. It works and it bridges, but
eight fishery fleets and five survey fleets that are really ages is not a pattern
to recommend, and the pin depends on `nonpar` selectivity being mean
standardized, which is an implementation detail rather than a contract.

## What already exists

Most of it. The model already predicts both quantities, per fleet, and reports
them:

    CAA      [pop, region, year, season, age, sex, fish_fleet]
    SrvIAA   [pop, region, year, season, age, sex, srv_fleet]

So the population dynamics, selectivity, catchability and mortality all stay
exactly as they are. What is missing is only the ability to *observe* those
arrays. This is an observation-layer feature, not a dynamics one.

Age-specific catchability needs no new routines either. A free catchability per
age is the same statement as one catchability times a non-parametric selectivity
at age, and `SPoRC` already estimates both. The aggregation is what removes the
age resolution, not the parameterization.

## Follow SAM's key matrices

SAM couples parameters through integer key matrices, `[n_fleets x n_ages]`, in
which equal entries share a parameter and `-1` excludes one. `keyVarObs` does
this for observation variances, `keyLogFpar` for survey catchability. One
mechanism covers every case that was going to need a separate spec string here:

| Intent | Key row |
|---|---|
| One variance per age | `1 2 3 4 5` |
| Variance by age group, as `smsR` does | `1 1 2 2 2` |
| One variance for the fleet | `1 1 1 1 1` |
| Age not observed by this fleet | `NA` |

`SPoRC` already thinks this way. `build_shared_spec_map()` and the `factor()`
maps throughout the setup functions are the same idea, and SAM's `-1` is
`SPoRC`'s `NA`. So the key matrices should be `SPoRC` maps, not a new
convention, and the existing identifiability check in
`check_spec_map_identifiable()` applies unchanged: a variance informed by fewer
than two observations is unbounded, and that is as true per age as it is per
year.

Concentrating a variance out analytically, which `smsR` does for catch through
`estSD`, is not something SAM offers. It costs no parameter and is the maximum
likelihood value given the residuals, so it is worth having as a fourth option
alongside the key matrix rather than as part of it.

## Proposed interface

### Catch at age

A fleet has either aggregated catch with compositions, or catch at age. Not
both: they are the same information stated twice, and fitting both double counts
it. This should be an error at setup, not a warning.

    Setup_Mod_Catch_and_F(
      ...,
      ObsCatchAA,          # [n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets]
      UseCatchAA,          # same shape, 0/1
      catchAA_units,       # "abd" or "biom", per fleet
      sigmaCAA_key,        # [n_ages, n_sexes, n_fish_fleets] integer map, NA excludes
      sigmaCAA_spec        # "fix", "est", or "concentrated"
    )

`catch_units` already exists and is reused. `UseCatchAA` governs which cells
are fished, replacing the `UseCatch` role, and the closure logic in
`model_population_dynamics.R:104` needs to read whichever data source the fleet uses.

### Index at age

    Setup_Mod_SrvIdx_and_Comps(
      ...,
      ObsSrvIdxAA,         # [n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets]
      UseSrvIdxAA,
      srv_q_key,           # [n_srv_fleets, n_ages], SAM's keyLogFpar
      sigmaSrvIdxAA_key,
      sigmaSrvIdxAA_spec
    )

with the fishery index mirroring it. When `srv_q_key` is supplied, catchability
is resolved per age from it and the fleet's selectivity is not estimated, since
the two are the same parameter twice over. That has to be enforced, or the model
is unidentified in a way nothing downstream will report.

### Correlation across ages

SAM's `obsCorStruct` is per fleet, one of independent, AR(1) across ages, or
unstructured. `SPoRC` already has `SrvIdx_LikeType` with an `"mvn"` option
and a supplied covariance, so the vocabulary exists; what is missing is an
*estimated* AR(1) across ages rather than a fixed matrix. Worth doing, but it is
separable from the rest and should come second. Independent first.

SAM's `fracMixObs`, a t(3) mixture for robustness, is a further increment and
should not hold up the first pass.

### Effort-driven fishing mortality

`smsR` writes fishing mortality as a seasonal pattern times a blocked age pattern
times an observed effort series, which is what makes an age-disaggregated catch
likelihood identifiable at all. Without it, four ages by thirty-nine years by two
seasons of free fishing mortality is exactly the number of observations, and the
model saturates.

    Setup_Mod_Catch_and_F(..., effort, use_effort, estimate_creep)

`Fmort` is built in exactly one place, `model_population_dynamics.R:110`, and
everything downstream consumes `Fmort` rather than its parameterization, so this
is one branch plus setup plumbing:

    Fmort[r,y,seas,f] = effort[r,y,seas,f] * exp(ln_F_mean[r,seas,f])

Two things to get right. `smsR` normalizes effort internally, so the raw series
and the fitted `Fagein` are on different scales; an estimated mean absorbs the
constant but anything seeding parameters from a reference fit will not. And
reference points need a decision: `Fmsy` under an effort-driven fishing mortality
means something different, and `refpts_*.R` currently assumes `Fmort` is free.

Note this is expressible today by fixing `ln_F_devs` at `log(effort)` and mapping
them off, which is how the sandeel case study does it. That is a legitimate way
to validate the feature before committing to the API.

## What this does not change

Population dynamics, selectivity, catchability, movement, tagging, growth,
recruitment. All of it is untouched. The feature reads arrays the model already
computes and compares them to data it currently has no way to accept.

## Work

Ordered so that each step is testable on the sandeel bridge before the next.

1. `combine`/`build` helpers for at-age observation error, next to
   `get_index_nLL` in `model_distributions.R`, following the `combine_idx_sd`
   and `build_idx_sd` pattern.
2. Catch at age: data arrays, validation, key mapping, likelihood block, the
   exclusivity error against compositions.
3. Index at age, both fishery and survey, with `srv_q_key` and the selectivity
   exclusivity error.
4. Concentrated variance as a third `spec` value.
5. Effort-driven fishing mortality.
6. AR(1) across ages, if wanted.

Every step holds the usual: roxygen on each argument with the allowed values
described, a section in `t_model_options.Rmd`, tests, simulation routing through
`sim_observations.R`, and `globals.R` for any new name.

## Tests worth writing

- An at-age likelihood on a fleet whose ages are all observed reproduces the
  aggregate likelihood when the key couples every age to one variance and the
  compositions are switched off. This is the sanity check that the two
  formulations agree where they should.
- A key with a single observation on some age errors through
  `check_spec_map_identifiable()`.
- Supplying both compositions and catch at age on one fleet errors.
- Supplying both `srv_q_key` and an estimated selectivity on one fleet errors.
- The concentrated variance equals the analytic maximum likelihood value of the
  residuals it is computed from.
- A structurally absent age, one not observed by a fleet, is excluded by `NA` in
  the key rather than relying on a likelihood tolerating a zero.
- The sandeel case study rebuilt on the real interface reproduces what the
  fleet-per-age formulation currently gets: spawning biomass to 6%, recruitment
  to 7%, the plus group to 16%, and the seeded comparison to 1.4%.

The last one is the real acceptance test. If the honest interface does not
reproduce what the workaround achieves, something in it is wrong.

## Open questions

- Should catch at age accept a total as well, so that a fleet can be fit at age
  while still reporting an aggregate for reference points? Reference points read
  catch, and if only the at-age data source is supplied they need to sum it.
- `smsR` applies fishing mortality to every age whether or not that age has a
  recorded catch, whereas a per-age observation naturally treats a missing cell
  as unfished. That difference is worth two cells on the sandeel bridge. The key
  matrix should distinguish "not observed" from "not fished", which one `NA`
  cannot do alone.
- Does an at-age fleet need its own `Wt_*` weight, and does weighting interact
  with a concentrated variance the way it does with an estimated one? It does,
  and the warning in `Setup_Mod_Weighting()` should extend to cover it.


## Addendum, 2026-08-27

The shapes above are as this note was first written. The data sources were since
generalized: every array has a sex dim, and each fleet names the dims
it sums over through `CatchAA_Type` and its counterparts, using the composition
vocabulary. Nothing promotes a short array; the declared shape is the only one
accepted. Reported standard errors, the choice of lognormal or normal, and where
the observation error comes from are all per fleet, matching what the aggregated
index data sources already allowed. The correlation options gained `"us"` and
`"2dar1"` alongside `"iid"` and `"1dar1"`, are chosen per fleet, and the
population-specific data sources have their own. See `vignette("t_model_options")`
and `vignette("ai_ices_style_assessments")` for the current interface.
