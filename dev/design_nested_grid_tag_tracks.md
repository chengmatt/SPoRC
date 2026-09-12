# A fine grid for tag tracks inside a coarse region assessment

Design note. Nothing here is implemented yet. Every number below was measured, and
the script that produced each is named where it matters.

## What this is for

Archival tags give daily positions at a few kilometres. An assessment runs on six
management areas. Collapsing each daily position into an area occupancy vector, which
is what the archived etag module does, throws away almost all of the movement
information: a fish that crosses West Yakutat in nine days and one that sits on a
single seamount for a year both read as "in area 5".

The proposal is one movement process on a fine grid, with the assessment's areas
defined as unions of fine cells. Tag tracks are fit at fine resolution. The population
dynamics stay at six areas, and their movement matrix is the fine transition matrix
aggregated over the areas.

## Why one process, not two tied together

`SPoRC` writes the movement rate out of region \eqn{j} into an adjacent region
\eqn{i} as

\deqn{q_{ij} = \theta_j / V_j + \max(\gamma_i - \gamma_j, 0)}

where

- \eqn{\theta_j} is the diffusion parameter for region \eqn{j}, units area per year,
  from `exp(2 * log_move_diffusion_pars)` through the diffusion design matrix
- \eqn{V_j} is the area of region \eqn{j}, `area_r[j]`, same area units
- \eqn{\gamma_i} is the preference (taxis) value of region \eqn{i}, dimensionless,
  from the preference design matrix
- the \eqn{\max} is the upwind form of `ctmc_diffusion_bounds`

On a square grid of side \eqn{h} this is \eqn{q = \theta / h^2} to each of four
neighbors, which is the finite volume discretization of two dimensional diffusion with
coefficient \eqn{\theta}. Mean squared displacement after time \eqn{t} is then
\eqn{4 \theta t} at **every** cell size, exactly. Checked at \eqn{h} = 3.131, 6.262,
12.524 and 25.048 km with \eqn{\theta} = 14 km²/day and \eqn{t} = 1 day: 56.0000 in
all four cases against \eqn{4 \theta t} = 56.

So \eqn{\theta} is a diffusion coefficient in real units, and the same value means the
same thing on the tag grid and on the management areas. One parameter serves both
resolutions, and the tag likelihood and the population dynamics are fit to the same
number rather than to two numbers that happen to be similar.

That also settles which grid the process belongs on. The six sablefish longline areas
are not squares, so \eqn{\theta_j / V_j} on them has no diffusion interpretation at
all; it is
whatever exchange rate makes the areas mix. The fine cells are squares, so on them the
expression means what it says. Put the process on the fine grid and derive the areas
from it, not the other way round.

The size of the difference, from `dev/scratch/nested_grid_lumping.R`: a 6 x 6 grid of
unit cells, four coarse regions of 3 x 3 cells each, \eqn{\theta} = 2, a preference
gradient that varies within each coarse region. Building the generator directly on the
four coarse regions with the same \eqn{\theta} and each region's mean preference gives
retention (the diagonal of the seasonal movement matrix) of 0.689 for region 1.
Aggregating the fine chain gives 0.576. The coarse only version holds fish 11 points
too long, because it treats a region as well mixed, and a fish that starts near a
boundary leaves far faster than the region average.

There is a second reason, already on file from the sandeel work: with no taxis,
equilibrium density is proportional to \eqn{V/\theta}, so a region varying diffusion
formula and a region varying preference formula are confounded in the equilibrium
distribution and separate only through the transient. A tag track is the transient,
sampled daily. This is the data source that breaks that confounding, and it cannot do
it at six areas.

## The aggregation step, and the one choice it forces

Write \eqn{P^{f}(\Delta)} for the fine transition matrix over one season,
\eqn{\exp(Q^{f} \Delta)}, and let region \eqn{R} be a set of fine cells. The coarse
movement matrix is

\deqn{P^{c}_{RS} = \sum_{c \in R} w_c \sum_{d \in S} P^{f}_{cd}, \qquad \sum_{c \in R} w_c = 1}

where

- \eqn{c}, \eqn{d} index fine cells, \eqn{R}, \eqn{S} index management areas
- \eqn{w_c} is the share of the region's fish sitting in cell \eqn{c} at the start of
  the season, dimensionless, summing to one within each region

The aggregate is free of \eqn{w} only when \eqn{\sum_{d \in S} P^{f}_{cd}} is the same
for every cell \eqn{c} in \eqn{R}, which is the lumpability condition for Markov
chains and does not hold for any real coastline. So \eqn{w} has to be chosen, and the
choice is the whole modeling content of the coupling. Three candidates:

**Area weights**, \eqn{w_c \propto V_c}. Says a region is uniformly occupied. This is
the well mixed assumption again, in weaker form, and it errs in the same direction as
the coarse only model.

**Equilibrium weights**, \eqn{w_c \propto \pi_c} for \eqn{c \in R}, with \eqn{\pi} the
stationary distribution of \eqn{Q^{f}}. This assumes only that the *shape within* a
region has settled, not that regions are in balance with each other, which is much
weaker. In the 6 x 6 test the equilibrium weighted aggregate reproduces the fine
chain's aggregated stationary distribution exactly (0.2555, 0.2445, 0.2445, 0.2555 by
both routes) while area weights give 0.2590 and 0.2410.

With no taxis, \eqn{\pi} has a closed form and needs no solve:

\deqn{\pi_c = (V_c / \theta_c) \big/ \sum_{c'} (V_{c'} / \theta_{c'})}

Verified to 5 decimals on a ring of 8 cells with random \eqn{\theta} and random areas.
With taxis there is no closed form; \eqn{\pi} is the null vector of \eqn{Q^{f}}, which
on a thousand cells is a sparse solve per gradient evaluation. Two ways out: run the
first version with pure diffusion at fine scale and preference only at the coarse
level, or fix \eqn{\pi} from a pilot fit and update it between phases.

**Abundance weights**, \eqn{w_c} = the model's own numbers at age in cell \eqn{c}.
Exactly right and not available, because it needs a fine scale numbers at age state.
Rejected.

## Sizing the fine grid

The six sablefish longline areas total 3,157,218 km²: Bering Sea 999,006, Aleutians
1,000,106, WGOA 329,822, CGOA 458,306, WY 189,088, EY/SE 180,890.

| Cell side | Cells over the full domain |
|---|---|
| 3.131 km (tag native) | 322,061 |
| 10 km | 31,572 |
| 25 km | 5,052 |
| 50 km | 1,263 |

Masking to sablefish depths should remove most of the Bering shelf and the deep
basin, which is where two thirds of that area sits. A 25 km grid masked to slope
depths lands in the low thousands, a 50 km grid in the low hundreds. Target 500 to
1500 cells.

Pick the cell side as an exact multiple of the tag grid so the daily probability
arrays aggregate by summation with no resampling: 8 x 3130.959 = 25,047.7 m,
16 x 3130.959 = 50,095.3 m. Each tag has its own grid origin, so the tags do not share
a lattice with each other; define the analysis grid once, then take zonal sums per
tag against the extent rebuilt from that tag's day one.

## What the tag files can supply

`hmm_results_<tag>_run<n>.RData` holds `probability.array`, the daily posterior
utilization distribution at 3.131 km, roughly 150 to 200 non-zero cells per day out of
64,512 in the crop.

Those posteriors are not observations. They come from a forward backward pass that
already assumed a movement model, and `model.parameters` states it: `D1` = 14, a
diffusion coefficient of 14 km²/day, with a numerical sub-division count of 2.856.
Fitting `SPoRC`'s \eqn{\theta} to them pulls it toward 14 km²/day whatever the light
and depth records say, because the smoothing put it there. Two tags smoothed with
different `D1` would disagree for no reason in the data.

The observation the fine likelihood wants is the daily emission probability, position
from light and depth against bathymetry, before any smoothing. `fish.data` in the same
file holds the inputs it is built from: `MinDepth` and `MaxDepth` with accuracy codes,
and `GeoLong` and `GeoLat` from the light record. Rebuilding those surfaces on the
analysis grid is a preprocessing job, not a modeling one, and it is the only route
that lets the assessment estimate diffusion rather than inherit it.

If rebuilding is too much for a first pass, use the Viterbi path as a hard position and
thin to one position per month. Losing information is safer than importing another
model's movement prior. Note that 14 km²/day is 5,110 km²/yr, so a starting value on
`SPoRC`'s scale is `log_move_diffusion_pars = 0.5 * log(5110) = 4.269`.

## Cost

Measured with `RTMB` on a 25 x 40 grid, 1000 cells and 3870 edges, a forward
recursion over 8 seasons at liberty with 4 sub-steps per season, one occupancy vector
per tag per season:

| Tracks | Sparse products | Tape build | One gradient |
|---|---|---|---|
| 10 | 320 | 0.7 s | 0.02 s |
| 50 | 1,600 | 3.8 s | 0.12 s |
| 150 | 4,800 | 11.1 s | 0.37 s |

`dev/scratch/nested_grid_tag_filter_cost.R`.

Linear in the number of tracks, and cheap enough to sit inside the objective. What
makes it cheap is never forming a transition matrix. Keep the generator as a vector of
edge rates and advance the occupancy vector \eqn{\alpha} by scatter:

```r
fl <- alpha[from] * q * dt                                   # flow along each edge
alpha <- alpha + as.vector(Sin %*% fl) - as.vector(Sout %*% fl)
```

with `Sin` and `Sout` fixed sparse 0/1 matrices, `[n_cells x n_edges]`, built once
outside the objective. `dt` must satisfy `max(outrate) * dt < 1` for the step to stay
positive; four sub-steps a season is ample at these rates.

Three things that do not work at this size:

- `Matrix::expm` on a thousand cells under AD. The dense result alone is 8 MB per
  season per age class.
- `mat_exp(..., expm_nsub)`, which is a dense `solve` followed by repeated squaring.
  Separately, that squaring loop only lands on the requested power when `expm_nsub` is
  a power of two: `nsub = 3` returns the 4th power of \eqn{(I - Q/3)^{-1}} and
  `nsub = 5` the 8th, off by 6.1e-2 and 1.0e-1 on a 3 region test. Worth fixing on its
  own account.
- Assigning into a large AD array inside the season loop. Keep \eqn{\alpha} a plain
  vector and use whole vector operations.

## What in SPoRC can be reused

The archived etag module in `SPoRC_etag_archive/` already does the hard part of the
tag likelihood: a per-tag ordered plan of move and observe steps that respects
`move_timing`, soft occupancy vectors, conditioning on a fishery determined endpoint so
movement is not pulled toward high F areas, survival weighting, terminal fate, and
mixing over unknown release population, age and sex weighted by abundance and release
platform selectivity. All of that is written against `n_regions` and is unchanged in
substance at `n_cells`. Two edits: the propagation step becomes the edge scatter above,
and `etag_obs_q` becomes `[n_etags, max_obs, n_cells]` sparse rather than dense.

What cannot be reused is `Get_Movement` itself. It allocates `Movement` and `Mrate` as
`[pop, from, to, year, seas, age, sex]`, which at 1000 cells is 7.2e9 entries, 58 GB.
The fine generator needs its own builder returning edge rates only, of length
`n_edges`, with year, age and sex structure applied as multiplicative offsets rather
than as extra array dims. For the same reason `ctmc_move_dat`, one row per population,
region, year, season, age and sex, becomes 7.2 million rows at 1000 cells; the fine
covariates have to be a lookup on cells alone, or cells by season.

Attribution and survival weighting stay coarse. Only the position state goes fine, and
it reads coarse numbers at age and coarse total mortality by looking up which area each
fine cell belongs to.

## Three routes, in order of commitment

1. **Plug in.** Fit the fine model outside `SPoRC`, aggregate to
   `[6, 6, season, age]`, pass it as `Fixed_Movement` with `use_fixed_movement = 1`.
   No package changes at all. Movement uncertainty does not reach the assessment, and
   the tags cannot respond to what the assessment learns about abundance.
2. **Prior.** Same fine fit, then put its \eqn{\hat\theta}, \eqn{\hat\gamma} and their
   covariance on the coarse CTMC parameters as a multivariate normal prior. Movement
   uncertainty propagates. The coarse generator is still built on polygons, so the
   over-retention above is untouched.
3. **Joint.** Fine generator inside the objective, coarse `Movement` by equilibrium
   weighted aggregation, tracks fit on the fine grid. This is the version that answers
   the question as asked, and the only one where the tags inform diffusion and
   preference separately.

Route 1 is a week and is worth doing first regardless, because it produces the
aggregation code and the analysis grid that route 3 needs, and it says immediately
whether the aggregated matrix looks anything like the movement estimates the
assessment produces on its own.

## Open questions

- The fine grid must be connected after the depth and land mask. A detached seamount
  group makes the chain reducible and \eqn{\pi} undefined per component. Check before
  anything else.
- Whether the equilibrium weights are refreshed inside the optimization or fixed
  between phases. Refreshing is a sparse null vector solve per gradient once taxis is
  on; fixing them makes the coarse matrix slightly inconsistent with the fine one.
- Whether age varying movement is worth keeping at fine scale. Sablefish age
  structure in movement is real, but 1000 cells by 30 ages is 30 times the tape, and
  the tag sample will not support age varying diffusion. Two or three age groups at
  fine scale, mapped up to the full age range at coarse scale, is the likely answer.
- Seasonal preference at fine scale needs a covariate that varies within a management
  area and within a year. Depth and slope are static; a spawning aggregation term would
  have to come from the fishery or the survey.

## What route 1 produced

Built and run on the Pacific cod archival tags, 2026-09-10. Code and results live with that
project, not here; this section records only what the run settles about the design.

**Diffusion transfers between grids, and that is what makes the aggregation work.** The fit
ran on a 1367 cell grid masked to bathymetry, and the aggregation on a 3432 cell grid
covering the three regions in full, at the same fitted diffusion. Cell size is identical, so
the local rate is identical, and only which cells exist changes.

**The weighting choice was a no-op on this dataset.** Equilibrium and area weights agreed to
3e-14. Equal area cells with no taxis make the fine chain's equilibrium uniform, so the two
are the same calculation. They separate only with a preference term or a cell-varying
diffusion, so the choice matters for the design and not for every application of it.

**A coarse generator over-retains, as predicted, and by about the size the toy suggested.**
Season retention for the western Gulf was 0.8206 from the aggregated fine chain against
0.8777 from a three region generator at the same diffusion, 5.7 points; the central Gulf
and the Bering were within 0.3 points, since both are large enough to be near well mixed at
this diffusion. The gap grows with the ratio of a season's displacement to the region's
width, so it is the small regions that need the fine grid.

**Cost was not a factor.** 11 tags, 1280 observations, 2464 daily steps, four sub-steps a
day: about 90 seconds to build the tape and fit, on the edge-scatter propagation. The
aggregation itself is 0.3 seconds, because only one propagation per region per season is
needed, never the fine transition matrix.

**The explicit step needs a sub-step count set from the rate, not chosen once.** A cell
losing more than its own occupancy in a step goes negative and the recursion diverges,
which arrives as an NaN objective partway through a search rather than as an error. The
fit also needs the diffusion bounded above for the same reason.

**Two things about the data turned out to matter more than anything in the method.** The
smoothed daily distributions in a geolocation output are not observations, because the
smoothing already assumed a movement model and records it. And the bathymetry the tag
geolocation was run on bounds where the fine grid can go, so it has to cover the assessment
regions before any of this is worth doing.

## Two corrections the Pacific cod run forced

**Equilibrium weights need a generator that settles somewhere, and a directional term does
not.** A preference surface linear in position has no restoring force, so its stationary
distribution is the downstream corner of the domain. With a fitted seasonal drift of about 4
km/day the equilibrium put 67% of the fish in one region in one season and 63% in another in
the next, and the aggregated matrix built on those weights was unusable, showing 65% of the
Bering moving to the western Gulf in a quarter. Area weights gave a sensible matrix from the
same generator. The rule is: equilibrium weights when the movement process has a stationary
distribution worth the name, area weights as soon as advection is on.

**Sub-stepping the propagation is not good enough, and the error is not where it looks.**
Advancing the occupancy vector by explicit steps of \eqn{(I + Q\,dt)} is first order, and on
this problem it was still 40% out on the faster season's diffusion at 32 steps a day, while
the slower season looked converged at 3%. Uniformization removes it: shift to
\eqn{P = I + Q/\lambda} with \eqn{\lambda} above the fastest exit rate, so \eqn{P} is an
ordinary transition matrix, and average \eqn{P^k} over a Poisson number of ticks. The term
count comes from the Poisson tail and is data, so it is a constant on the tape. Agreement
with a dense matrix exponential was 1e-13, with or without mortality on the diagonal.

The three endpoint cases are worth stating once. A scheduled pop-up is a movement
observation and the tag is conditioned on survival. A fleet recovery is where the fleet was
fishing, so it goes into both the numerator and the denominator, restricting which paths
count without being evaluated as movement, and the fishing mortality at that cell cancels
between the two, so the conditional form needs no F. The joint form keeps the capture in the
numerator only and lets the tag inform F, at the price of not being able to count the same
recapture in a conventional tag likelihood.

## Which scale the rates should be defined at

Diffusion transfers between grids and, in the rate form of preference SPoRC had until
2026-09-11, preference did not. The measurements below are of that form; the section after
this one records the form that makes both transfer.

\eqn{\theta} is a diffusion coefficient: \eqn{\theta/V = D/h^2} on a square cell, so mean
squared displacement is \eqn{4\theta t} at every cell size, verified exact. A preference
coefficient is a raw rate difference along an edge, and the covariate gradient between
neighbours grows with how far apart they are, so the drift it produces goes as \eqn{h^2}.
Measured on a chain covering a fixed 800 km with the same coefficient of 0.002 per km:

| Cell side | Net drift |
|---|---|
| 12.5 km | 0.31 km/day |
| 25 km | 1.25 km/day |
| 50 km | 4.91 km/day |
| 100 km | 12.12 km/day |

Rescaling the coefficient as \eqn{v/h^2} holds the drift at 1.25 km/day at all four sizes.
So a preference coefficient estimated at one resolution means nothing at another, while a
diffusion coefficient means the same thing everywhere.

Three consequences.

**Aggregate the rates, not the covariate.** Averaging a covariate over a region and then
applying the link is not the same as applying the link cell by cell and aggregating. On the
cod grid the Bering's mean depth of 70 m hides a range of 11 to 358 m and the western Gulf's
106 m hides 36 to 321, which is exactly the within-region structure that drives movement.

**A coarse covariate is a special case of a fine one.** A region level or region by season
covariate, a heatwave flag or a survey index, is broadcast to the cells of that region with
no approximation. So nothing usable at the coarse scale is lost by defining finely, while
anything varying within a region is unavailable the other way.

**The conventional tag component is not lost.** It reads the movement and rate arrays, which
is what the aggregation produces, so nothing about it changes. There is a real point next to
it though: a conventional release is a position, not a region average, and the release
apportionment currently spreads a cohort over its whole region by abundance and selectivity.
With a fine grid each cohort can start in its own cell, which is both closer to the data and
already in it.

## Preference made scale-free, 2026-09-11

The rate form `q = theta_j / V_j + max(gamma_i - gamma_j, 0)` divides diffusion by area and
not the contrast, and that alone is why the drift went as \eqn{h^2}. Writing the contrast as
a multiplier of the edge's diffusive rate,

    q_ji = theta_j / V_j * (1 + max(gamma_i - gamma_j, 0)),

makes gamma dimensionless. On a 1-D lattice with a linear potential, the slope of log
equilibrium density per unit distance, whose continuum value is 0.02:

| Form | 1 km cells | 8 km cells |
|---|---|---|
| rate form, same gamma | 0.0198 | 0.3098 |
| contrast form, same gamma | 0.0198 | 0.0186 |
| exponential, theta exp(d / 2) | 0.0200 | 0.0200 |

The contrast form is now what SPoRC's generator does (`R/model_movement.R`) and what every
fine grid script does. The exponential form is exact for the equilibrium at any cell size
because it satisfies detailed balance, and would be a fourth bound option if ever wanted. A
step contrast across a region boundary is also nearly scale-free on the fine grid: the
seasonal exchange it produces moves from 0.0773 at 1 km cells to 0.0788 at 16 km.

What the form does not fix, and nothing can, is a coarse cell many diffusion lengths wide.
With theta of 5 km2/day and a 180-day season the diffusion length is 60 km. A 256 km region
lumped from the fine chain exchanges 6.6% per season under diffusion alone, where a two-state
chain handed the same theta gives 1.4%. The coarse rate that reproduces the fine exchange is
5.2 times the copied one, and it grows as the square root of theta (log slope 0.56) rather
than linearly. A step contrast of 2 on the fine grid needs a coarse contrast of 8.6 and a
coarse theta of 21.6 to reproduce the season; a step of 20 needs 35.2 and 8.8. So coarse
rates are effective rates set by theta, region size and season length together. They come
from lumping the fine seasonal matrix, never from copying parameters. That is why route 3
recomputes the coarse matrix from the fine generator, and why a movement covariate has to
enter the fine generator in every year, with conventional-tag years reading its lumped
consequence, rather than living in a coarse model for some years and a fine one for others.

A uniformization defect was found the same day. The shift was sized for diffusion alone, and
the region-preference fit of 2026-09-10 put contrasts of 12 to 35 times the diffusive rate on
boundary edges, so \eqn{P = I + Q/\lambda} had entries below -1 there and the series diverged:
mass of -9.8e9 after 30 days from the busiest cell. The shift now sits on the tape as 1.25
times the busiest exit rate, floored at a design bound that fixes the chunk and term counts,
and every fit reports `shift_ratio`. Results from the rate-form region-preference fit are
superseded by the refit recorded in the fine grid README.

## Decision, 2026-09-11: the coarse-first construction is the method

The assessment stays on its three regions with the coarse generator Q as the model of movement.
The tracks read Q spread over the fine cells, every cell of region S jumping into region R at
rate Q[S, R] and landing evenly over R's cells, plus diffusion within regions at one rate per
season. That fine generator lumps back to Q exactly at any interval and from any start (Z M = Q Z),
so no calibration or penalty is needed and a covariate on movement lives on Q alone, whatever
data the years it varies in hold. The price is a region without memory: a fish that just crossed
is as likely to cross back as a long-time resident, and a crossing lands anywhere. Under a truth
with directed movement at the boundary the tracks-only self-test on the cod grid puts the error
on a seasonal fraction at 0.25 to 0.50; that is the recorded caveat, not a defect to fix.

Everything else tried between 2026-09-10 and 2026-09-11 is archived under
`Pcod_Spatial/R/fine_grid/explored/` with its results: the standalone track fits and their
handoff (rejected by the assessment data), the fine generator inside the objective with
covariates (route 3; misspecified when movement is directed at the boundary), local crossing
calibrated on the tape (route 4b; a diffusion ceiling and solver gradients), local crossing
with estimated permeabilities and a penalty on the season matrices (route 4c; works, but
agreement holds only at the season step from an even start), an advective-only boundary
(route 4d; fails when diffusion crosses), and zones (exact, but more coarse states). The main
line is `01_fine_grid.R` to `05_simulate_tracks.R` plus the drivers, described in that folder's
README. The SPoRC package itself is unchanged by this beyond the scale-free preference form
recorded above.

## Reverted, 2026-09-12: the preference form

The contrast form of preference (theta / area times one plus the contrast) was reverted at
Matt's request; the package is back to the rate form it had before 2026-09-11 and the
scale-invariance test was removed. The coarse-first construction does not need it: the fine
generator is built from whatever coarse generator the assessment produces and lumps back to it
exactly under either form. The contrast form served only the fine-first idea, which was set
aside. Under the cod configuration the two forms are different models: objective 3539.4 against
3546.5, terminal spawning biomass 102.2 against 105.2 kt, up to 0.11 on a seasonal crossing
fraction. The fit made under the contrast form is kept as
`output/fine_grid/cod_move_estimated_contrastform.RDS` beside the rate-form comparator.
