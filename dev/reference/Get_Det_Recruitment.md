# Deterministic Recruitment

Computes deterministic recruitment by population and region using either
a mean recruitment model or a Beverton-Holt stock-recruitment
relationship.

## Usage

``` r
Get_Det_Recruitment(
  recruitment_model,
  rec_dd,
  y,
  rec_lag,
  R0,
  rec_region_prop,
  rec_seas_prop,
  h,
  n_pop,
  n_regions,
  n_ages,
  n_fish_fleets,
  WAA,
  MatAA,
  natmort,
  SSB_vals,
  Movement,
  sgl_seas_spawning_movement,
  stray_rate,
  do_recruits_move,
  t_spawn,
  init_F,
  dmr,
  fish_sel,
  ret_sel,
  n_seas,
  spawn_seas,
  natal_region,
  seasdur,
  sexratio_f,
  Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0
)
```

## Arguments

- recruitment_model:

  Integer flag specifying the recruitment model:

  - `0` Mean recruitment

  - `1` Beverton-Holt recruitment with steepness

  - `2` Ricker recruitment with steepness

- rec_dd:

  Integer flag specifying the density dependence structure:

  - `0` Local density dependence (population or region specific)

  - `1` Global density dependence (shared across regions; only valid
    when `n_pop = 1`)

- y:

  Current model year index.

- rec_lag:

  Recruitment lag (in seasons) between spawning and recruitment. `1` is
  the classic lagged case: recruitment uses `SSB_vals` from `rec_lag`
  seasons prior. `0` is age-0 recruitment: recruitment uses the SAME
  year's SSB (`SSB_vals[,,y]`). The caller is responsible for supplying
  that value already computed from survivors only (i.e. before this
  year's recruits exist) when `rec_lag = 0`, see `SPoRC_rtmb.R`,
  `Simulate_Population.R`, and `Do_Population_Projection.R` for how each
  population-dynamics loop does this.

- R0:

  Numeric vector (`n_pop`) of unfished recruitment by population.

- rec_region_prop:

  Matrix (`n_pop × n_regions`) giving the proportion of recruitment
  allocated to each region.

- rec_seas_prop:

  Matrix (`n_pop × n_seas`) giving seasonal recruitment proportions.
  When `rec_lag = 0`, must be zero for every season before `spawn_seas`
  (age-0 recruits can't predate the spawning event that produced them),
  validated at setup by `Setup_Mod_Rec`/`Setup_Sim_Rec`.

- h:

  Matrix (`n_pop × n_regions`) of Beverton-Holt steepness values.

- n_pop:

  Number of populations.

- n_regions:

  Number of spatial regions.

- n_ages:

  Number of age classes (including the plus group).

- n_fish_fleets:

  Integer. Number of fishery fleets.

- WAA:

  Array (`n_pop × n_regions × n_seas × n_ages`) of weight-at-age.

- MatAA:

  Array (`n_pop × n_regions × n_seas × n_ages`) of maturity-at-age.

- natmort:

  Array (`n_pop × n_regions × n_ages`) of natural mortality.

- SSB_vals:

  Array (`n_pop × n_regions × n_years`) of spawning biomass.

- Movement:

  Array (`n_pop × origin × destination × n_seas × n_ages`) giving
  seasonal movement probabilities.

- sgl_seas_spawning_movement:

  Array (`n_pop × origin × destination × n_ages`) describing spawning
  movement when a single season is used and `n_pop > 1`.

- stray_rate:

  Numeric vector of stray rates by population.

- do_recruits_move:

  Indicator for whether recruits move in their first year.

- t_spawn:

  Fraction of the spawning season that occurs before spawning.

- init_F:

  Array (`n_regions × n_seas × n_fish_fleets`) of initial fishing
  mortality.

- dmr:

  Array (`n_regions × n_seas × n_fish_fleets`) of initial (first year)
  discard mortality.

- fish_sel:

  Array (`n_pop x n_regions × n_seas x n_ages x n_fish_fleets`) of total
  fishery selectivity.

- ret_sel:

  Array (`n_pop x n_regions × n_seas x n_ages x n_fish_fleets`) of
  retained fishery selectivity.

- n_seas:

  Number of seasons per year.

- spawn_seas:

  Season index in which spawning occurs.

- natal_region:

  Integer vector (`n_pop`) mapping each population to its natal region.

- seasdur:

  Numeric vector (`n_seas`) giving seasonal durations as fractions of a
  year.

- sexratio_f:

  Matrix (`n_pop × n_regions`) giving female recruitment proportions.

## Details

Recruitment is distributed spatially using regional recruitment
proportions and seasonal recruitment timing. When Beverton-Holt
recruitment is used, unfished spawning biomass per recruit (\\S_0\\) is
calculated internally by projecting a single recruit through the full
seasonal population dynamics, including movement and mortality.

Two recruitment formulations are supported.

\*\*Mean recruitment\*\*

When `recruitment_model = 0`, recruitment is constant:

\$\$R\_{p,r} = R\_{0,p} \times RecProp\_{p,r}\$\$

where recruitment is distributed spatially according to
`rec_region_prop`.

\*\*Beverton-Holt recruitment\*\*

When `recruitment_model = 1`, recruitment follows the Beverton-Holt
relationship:

\$\$ R = \frac{4hR_0SSB}{(1-h)S_0 + (5h-1)SSB} \$\$

where:

- \\SSB\\ is spawning biomass lagged by `rec_lag` seasons (or, when
  `rec_lag = 0`, the current year's own spawning biomass, see the
  `rec_lag` parameter above)

- \\S_0\\ is unfished spawning biomass per recruit

- \\h\\ is steepness

\*\*Ricker recruitment\*\*

When `recruitment_model = 2`, recruitment follows the Ricker
relationship, written in the depletion form used by the EBS pollock
assessment (2024, `SrType = 1`):

\$\$ R = R_0 \frac{SSB}{S_0} \exp\left(\alpha \left(1 -
\frac{SSB}{S_0}\right)\right), \quad \alpha =
\log\left(\frac{4h}{1-h}\right) \$\$

The curve passes through \\(S_0, R_0)\\ by construction. Note that
\\\alpha\\ is set so the Ricker carries the same compensation ratio as a
Beverton-Holt at the same \\h\\, rather than by the textbook definition
\\R(0.2 S_0) = h R_0\\. Steepness is therefore not interchangeable
between the two curves: the Ricker here gives \\R(0.2S_0)/R_0 =
0.2(4h/(1-h))^{0.8}\\, which exceeds \\h\\ and is not bounded by 1.

\\S_0\\ (and the age-composition of spawning biomass per recruit more
generally) does not depend on `rec_lag`; it is a pure per-recruit,
equilibrium quantity. The recruit age class (the first age) is always
included in the sum; when `rec_lag = 0`, maturity at that age is
required to be exactly zero (validated at setup by
`Setup_Mod_Biologicals`/`Setup_Sim_Biologicals`), so it contributes
nothing regardless.

Spawning biomass per recruit (\\S_0\\) is computed internally by
projecting a single recruit through all ages and seasons under both
unfished and fished conditions. The algorithm:

1.  Allocates a recruit across regions and seasons.

2.  Applies seasonal movement.

3.  Applies natural and fishing mortality, where fishing mortality is
    decomposed into retained (\\F \cdot sel \cdot ret\\) and dead
    discard (\\F \cdot sel \cdot (1 - ret) \cdot dmr\\) components.

4.  Computes spawning biomass during the spawning season.

5.  Solves the plus group analytically using annual transition matrices.

When multiple populations are modeled, recruitment for each population
depends on spawning biomass in its natal region. Contributions from
other populations are scaled by the specified stray rates.
