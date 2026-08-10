# Overview of Model Options

SPoRC models are assembled through a pipeline of `Setup_Mod_*` functions
(for the estimation model) and `Setup_Sim_*` functions (for the
operating model / simulation). Each function appends to an `input_list`
(or `sim_list`) that accumulates data, parameter starting values, and
RTMB factor maps. This vignette catalogues every user-facing option
across that pipeline. Where both an integer code and a string alias
exist, either is accepted. Strings are generally preferred for
readability. When `n_sexes > 1`, the first sex index is always female.
When `n_pop > 1`, populations follow the order in which they are
defined. Throughout this document, the notation
`[p × r × y × τ × a × s]` is shorthand for array dimensions: population,
region, year, seasons, age, sex.

## Parameter sharing conventions

Many hyperparameters throughout the pipeline (process-error standard
deviations, catchability, natural mortality, tag reporting rates, …) are
controlled by a `_spec` string argument that follows the same
convention: a single named dimension (e.g. region, season, fleet, sex)
is either estimated independently for every level, shared across some
subset of levels, or fixed. Concretely, for a hyperparameter varying
across dimensions abbreviated `r` (region), `seas` (season), and `f`
(fleet), the accepted strings are:

| String | Meaning |
|----|----|
| `"est_all"` | Independent parameter for every combination of dimensions |
| `"est_shared_r"` | Shared across regions; independent per remaining dimension |
| `"est_shared_seas"` | Shared across seasons |
| `"est_shared_f"` | Shared across fleets |
| `"est_shared_r_seas"`, `"est_shared_r_f"`, `"est_shared_seas_f"` | Shared across the named pair of dimensions |
| `"est_shared_r_seas_f"` | A single parameter shared across everything |
| `"fix"` | Held at its starting value (not estimated) |

Which dimension abbreviations are valid, and in what order they can be
combined, depends on the specific parameter (see each section below for
its abbreviations, e.g. `sigmaC_spec` additionally has a `y` (year)
dimension). `sigmaF_spec`, `sigmaC_spec`, `sigmaR_spec`, and
`Fdev_rho_spec` all follow this convention via a single shared internal
helper (`build_pe_map`/`build_shared_spec_map`), so they behave
identically. Other `_spec` arguments using the same string vocabulary
(e.g. `sigmaC_pop_spec`, and selectivity’s `fish_fixed_sel_pars_spec`,
`fishsel_pe_pars_spec`, `fish_sel_devs_spec`, and their retention/survey
equivalents) are implemented separately with their own hand-written
logic, and some additionally support a fleet-sharing escape hatch,
`"est_shared_f_<n>"`, which copies the sharing structure of fleet `<n>`
wholesale rather than collapsing a dimension.

------------------------------------------------------------------------

## Spatial and Demographic Structure

Everything starts with
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
(estimation model) or
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md)
(operating model; simulation). These functions define the skeleton that
every subsequent setup call validates against.

**Core dimensions**

| Argument | Type | Description |
|----|----|----|
| `years` | integer vector | Calendar years included in the assessment (length determines $`n_y`$) |
| `ages` | integer vector | Modelled age classes; the final element is the plus group |
| `lens` | numeric vector or `NULL` | Length-bin midpoints; `NULL` disables length-based features |
| `n_regions` | integer | Spatial regions |
| `n_sexes` | integer | `1` (aggregated) or `2` (sex-structured) |
| `n_fish_fleets` | integer | Fishery fleets |
| `n_srv_fleets` | integer | Survey fleets |

**Multi-population structure**

SPoRC decouples biological populations from spatial regions. Populations
have their own stock-recruit dynamics and natal assignment; regions
define the arena through which all populations move.

| Argument | Type | Description |
|----|----|----|
| `n_pop` | integer | Number of biological populations (default `1`) |
| `natal_region` | integer vector (length `n_pop`) | Maps each population to its home region. Inferred automatically when `n_pop == n_regions` (one-to-one) or `n_regions == 1` (all share region 1). Must be supplied explicitly otherwise |

When `n_pop > n_regions`, multiple populations share a natal region, a
contingent-population or mixed-stock structure. When
`n_pop < n_regions`, the model represents a single population spread
across more regions than there are spawning origins.

**Seasonal sub-stepping**

| Argument | Type | Description |
|----|----|----|
| `n_seas` | integer | Seasons per year (default `1`) |
| `seasdur` | numeric vector (length `n_seas`) | Duration of each season as a fraction of the year (must sum to 1). Defaults to equal-length seasons |

Within each year, processes execute in sequence: recruitment → movement
→ mortality (with age advancement at the end of the final season).

**Other dimension controls** (EM only)

| Argument | Description |
|----|----|
| `n_proj_yrs_devs` | How many projection-year slots to pre-allocate for deviation parameters (`ln_RecDevs`, `move_devs`, selectivity deviations). Default `0` |
| `store_config` | If `TRUE`, archives all `Setup_Mod_*` arguments inside `input_list$config` for reproducibility |

**Additional simulation dimensions** (OM only)

| Argument | Description |
|----|----|
| `n_sims` | Number of Monte Carlo replicates |
| `n_yrs` | Projection horizon (years) |
| `n_obs_ages` | Number of observed age bins (can differ from `n_ages` when compositions pool ages differently) |
| `run_feedback` | `TRUE` for closed-loop MSE; `FALSE` (default) for open-loop simulation |
| `feedback_start_yr` | First year of feedback when `run_feedback = TRUE` |

------------------------------------------------------------------------

## Population Initialisation

Controlled via
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
(argument `init_age_strc`) and
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md).

SPoRC offers four methods for deriving equilibrium numbers-at-age, plus
a fifth, non-equilibrium option. The equilibrium methods project a
constant recruitment ($`R_0`$, or a separate `ln_rinit` scalar when
`use_rinit = 1`) through seasonal mortality, fishing, and (optionally)
movement until the age structure stabilises. After the equilibrium is
computed, multiplicative log-scale initial deviations (`ln_InitDevs`)
are applied to ages $`2, \ldots, A`$.

When `use_rinit = 1`, the regional recruitment scalar feeding the
equilibrium calculation is bias-corrected:
$`\zeta_{p,r}\,\text{rinit}_p\exp(-\sigma_{\text{Rec,early}}^2/2)`$,
treating `ln_rinit` as the median of the assumed lognormal recruitment
process and converting to the corresponding mean before it is used as a
deterministic (non-stochastic) equilibrium seed. This mirrors the same
correction applied to
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)’s
equivalent `rinit_input` pathway during simulation, so fitted and
simulated equilibria remain on a consistent scale.

#### Equilibrium solution method (`init_age_strc`)

| Value | String | How the equilibrium is solved |
|----|----|----|
| 0 | `"iterative"` | Brute-force iteration: runs the full seasonal cycle $`n_\text{ages} \times 10`$ times |
| 1 | `"scalar_no_move"` | Closed-form geometric series assuming no movement at any age. Plus group: $`N_{A^+} = N_{A-1} e^{-Z_{A-1}} / (1 - e^{-Z_{A^+}})`$ |
| 2 | `"matrix"` | Builds seasonal transition matrices $`\mathbf{T}_a = \prod_\tau \mathbf{M}_{a,\tau} \mathbf{S}_{a,\tau}`$ that combine movement and survival, then solves $`\mathbf{N}_{A^+} = (\mathbf{I} - \mathbf{T}_{A^+})^{-1} \mathbf{T}_{A-1} \mathbf{N}_{A-1}`$. **Default** |
| 3 | `"scalar_plus_only"` | Hybrid: uses the matrix approach for ages below the plus group but switches to the scalar geometric-series for the plus group itself |
| 4 | `"free"` | No equilibrium at all: `ln_InitDevs` *are* the initial log numbers-at-age (ages 2+, apportioned by sex ratio), with age 1 still taken from recruitment |

**Recommendations.** Use the default `"matrix"` (2) for spatial models;
`"scalar_no_move"` (1) is fine, and marginally cheaper, for
single-region models. Use `"free"` (4) when the initial age structure
carries no information about $`R_0`$ and should not be pulled toward an
equilibrium, matching assessments in which initial numbers-at-age are
freely estimated parameters. Note two consequences of `"free"`: the
initial condition becomes independent of `init_F_par` and of `ln_rinit`,
and the deviations are on the scale of log-*numbers* rather than
log-ratios about an equilibrium, so any penalty applied via
`equil_init_age_strc` acts as a prior on log initial abundance. Pair
`"free"` with `equil_init_age_strc = "equil"` (0) if no such prior is
wanted, or with `InitDevs_pen_center = "own_mean"` (see Recruitment
below) to penalize only the roughness of the initial age structure
rather than its level.

#### Initial deviation structure (`equil_init_age_strc`)

| Value | String | Which ages receive initial deviations |
|----|----|----|
| 0 | `"equil"` | None, strict equilibrium |
| 1 | `"stoch_no_plus"` | All ages except the plus group. **Default** |
| 2 | `"stoch_all"` | All ages including the plus group |
| 3 | `"stoch_shared_ages"` | User-defined age sharing via `init_age_devs_shared` |

`init_F_par` (array `[n_regions × n_seas × n_fish_fleets]`) optionally
introduces fishing mortality into the equilibrium calculation, producing
a fished initial condition. `init_F_form` sets what it means: `"prop"`
treats it as a proportion of the mean fishing mortality (inverse-logit
scale, bounded to (0,1)), so the initial age structure moves with mean
F; `"abs"` treats it as an absolute rate (log scale) that is independent
of mean F. `init_F_spec` (`"fix"` or `"est"`) sets whether it is
estimated, independently of the form. Prefer `"abs"` when the historical
fishing mortality that shaped the initial condition is distinct from the
mean F of the modelled period; under `"prop"` one parameter both
depletes the initial age structure and scales the F series, and because
catch constrains only their product the optimizer can fit catch equally
well with a smaller, harder-fished stock. The older `init_F_prop`
argument is still accepted and is converted to the `"prop"` form.

------------------------------------------------------------------------

## Recruitment

Controlled via
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
(EM) or
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)
(OM).

#### Stock-recruit function (`rec_model`)

| String | Description |
|----|----|
| `"mean_rec"` | Estimate mean $`\ln R_0`$ with annual deviations; no SSB feedback |
| `"bh_rec"` | Beverton-Holt: $`R = 4hR_0 \cdot \text{SSB} / \left[(1-h)S_0 + (5h-1)\text{SSB}\right]`$. Unfished spawning biomass per recruit ($`S_0`$) is computed internally by projecting a single recruit through all ages and seasons with movement |
| `"ricker_rec"` | Ricker in depletion form: $`R = R_0 (S/S_0)\exp(\alpha(1 - S/S_0))`$ with $`\alpha = \log(4h/(1-h))`$ (the “Dorn form” used by the EBS pollock assessment). The curve passes through $`(S_0, R_0)`$ and shares the Beverton-Holt’s compensation ratio at a given $`h`$, not the textbook $`R(0.2S_0) = hR_0`$ definition, so steepness values are not interchangeable between `"bh_rec"` and `"ricker_rec"`, and a steepness prior calibrated for one should not be reused for the other |

**`SR_ref_yr`** sets the year index whose biological inputs (WAA,
maturity, $`M`$, movement) feed the unfished
spawning-biomass-per-recruit calculation, and hence $`S_0`$ and the
scale of the stock-recruit curve. Default `1` (first model year, the
long-standing behaviour); set `length(years)` to condition the curve on
terminal biologicals, a convention several assessments use. With
time-varying weight-at-age the two choices give different $`S_0`$; pick
whichever the assessment you are bridging or the reference-point
convention you follow uses, and keep it consistent with how reference
points are computed.

#### Density dependence scope (`rec_dd`)

| String | When to use |
|----|----|
| `"local"` | SSB and $`S_0`$ computed per population and/or per region. **Required** when `n_pop > 1` with BH recruitment |
| `"global"` | SSB summed across all regions before entering the SRR. Single-population only |

#### Spawning timing

| Argument | Description |
|----|----|
| `spawn_seas` | Season index in which spawning occurs |
| `t_spawn` | Fraction of the spawning season elapsed before spawning occurs (0 = start, 0.5 = midpoint) |
| `rec_lag` | Delay (in seasons) between spawning and recruitment entry. `1` (default): recruitment driven by SSB from `rec_lag` seasons prior, may enter in any season. `0`: age-0 recruitment, recruitment driven by that *same* year’s own SSB. Since that SSB isn’t known until `spawn_seas` is reached, recruits may only enter in `spawn_seas` itself or a later season in the same year (`rec_seas_prop` must be zero for every season before `spawn_seas`), and the recruit age class must have zero maturity everywhere |

#### Recruitment allocation

| Argument | Description |
|----|----|
| `rec_region_prop_spec` | How regional recruitment proportions are estimated/fixed |
| `rec_seas_prop_spec` | How seasonal recruitment proportions are estimated. Default `"fix"` (all recruitment enters in season 1) |
| `sexratio_spec` | Sex ratio at recruitment estimation. Default `"fix"` (equal or user-supplied) |
| `sexratio_blocks` | Time blocks for sex-ratio parameters, specified as `"none_Pop_<p>_Region_<r>"` or `"Block_<b>_Year_<s>-<e>_Pop_<p>_Region_<r>"` |

#### Recruitment and Initial age deviations

| Argument | Description |
|----|----|
| `RecDevs_spec` | Estimation structure for annual $`\ln(\text{RecDevs})`$ |
| `InitDevs_spec` | Estimation structure for initial age deviations |
| `sigmaR_spec` | Recruitment variability: `"est_all"`, `"fix_early_est_late"`, `"est_shared_all"`, `"fix"` |
| `sigmaR_switch` | Year index separating the early and late $`\sigma_R`$ periods |
| `do_rec_bias_ramp` | Methot & Taylor bias adjustment (0 = off, 1 = on) |
| `bias_year`, `max_bias_ramp_fct` | Bias ramp breakpoints and maximum correction factor |
| `dont_est_recdev_last` | Number of terminal-year recruitment deviations to fix at zero |
| `init_age_devs_shared` | Integer vector of length `n_ages - 1` specifying age-sharing for `ln_InitDevs`. Positions with the same value share a single estimated parameter (e.g. `c(1:42, rep(42, 9))`). Required when `equil_init_age_strc = 3`; `NULL` (default) uses standard behaviour |
| `use_rinit` | 0 = population initialised using `ln_global_R0` (default); 1 = separate `ln_rinit` used for initialisation, with `ln_global_R0` governing only the recruitment relationship. |

#### Deviation penalty centering (`RecDevs_pen_center`, `InitDevs_pen_center`, `Fdev_pen_center`)

The recruitment, initial age, and (iid) fishing mortality deviation
penalties can each be centred in one of two ways:

| String | Penalty mean | What it constrains |
|----|----|----|
| `"fixed"` | The asserted prior mean: zero, or the bias-corrected $`-\sigma_R^2/2`$ for recruitment under the bias ramp. **Default** | Both the level and the spread of the deviations |
| `"own_mean"` | The mean of the estimated deviations themselves | Only their spread; the level is left free (a sum of squares about the mean, matching assessments whose deviation vectors sum to zero) |

**Recommendations.** `"fixed"` is the statistically coherent choice when
$`\sigma_R`$ is estimated or the deviations are treated as random
effects: the penalty is then a genuine distributional assumption, and
the bias ramp machinery depends on the mean being asserted. `"own_mean"`
is primarily a *bridging* device: it reproduces assessments where the
mean parameter (`ln_global_R0`, `ln_F_mean`) carries the level and the
deviations carry only shape, so the level is not penalized twice. Two
cautions: under `"own_mean"` the deviations’ level must be pinned
elsewhere (an $`R_0`$ prior, a fixed deviation, or informative data) or
the likelihood is flat along it, and `RecDevs_pen_center = "own_mean"`
cannot be combined with `do_rec_bias_ramp = 1` (the $`-\sigma^2/2`$
offset is meaningless once the mean is estimated rather than asserted;
setup errors out).

#### Recruitment level penalty

Separate from the deviation penalty, an optional penalty on the log
recruitment *series itself*:

| Argument | Description |
|----|----|
| `Use_rec_level_pen` | 0 (default) / 1 toggle |
| `rec_level_pen_sigma` | Standard deviation of the penalty. A sum of squares with weight $`w`$ corresponds to $`\sigma = 1/\sqrt{2w}`$. Default 1 |
| `rec_level_pen_center` | `"own_mean"` (default; penalizes only the series’ variability) or `"fixed"` (centres on zero) |
| `rec_level_pen_yrs` | Calendar years the penalty applies over; `NULL` (default) = all years |

**Rationale.** Under a stock-recruit relationship the deviations are
residuals about the predicted curve, so a model that also wants the
*realized* recruitment series to stay regular has nowhere else to say
so. This penalty is that second, independent statement, and reproduces
the recruitment regularity penalties several existing assessments carry
(e.g., a penalty on log recruitment variability). Leave it off unless
recruitment in data-poor years is wandering unreasonably; it is a tuning
penalty, not a probability model, and it will shrink genuine recruitment
variability if over-weighted (hence the `rec_level_pen_sigma` is left as
data rather than an estimated parameter).

#### Recruitment penalty weighting (`Setup_Mod_Weighting()`)

| Argument | Description |
|----|----|
| `Wt_Rec` | Weight on the recruitment deviation penalty. Scalar (default 1), or an array `[n_pop × n_regions × n_est_rec_devs]` for per-deviation weighting; note the third dimension follows `ln_RecDevs` (moved by `dont_est_recdev_last` and `n_proj_yrs_devs`), not the number of years |
| `Wt_Init_Rec` | Weight on the initial age deviation penalty, `[n_pop × n_regions × (n_ages - 1)]` or scalar. `NULL` (default) inherits a scalar `Wt_Rec`; must be supplied explicitly when `Wt_Rec` is an array, since the two penalties are dimensioned differently |

A per-deviation weight of zero removes that deviation from the penalty
while it remains estimated, which is how a stock-recruit relationship is
fit over a chosen window of years while recruitment stays effectively
free elsewhere. That is distinct from `dont_est_recdev_last`, which
removes the deviations themselves (recruitment reverts to the
deterministic prediction), and from mapping a deviation off (fixed *and*
unpenalized; see [Which deviations are
penalized](#which-deviations-are-penalized)). Beware the $`\sigma_R`$
interaction: deviations excluded from the penalty no longer inform an
estimated `ln_sigmaR`, so weight-based windows are safest with
$`\sigma_R`$ fixed.

#### Steepness priors

When a stock-recruit curve is used (`rec_model = "bh_rec"` or
`"ricker_rec"`), steepness ($`h`$) can be penalised with a Beta
distribution scaled to \[0.2, 1\] by default:

| Argument | Description |
|----|----|
| `Use_h_prior` | 0 = no prior, 1 = apply prior |
| `h_prior` | Data frame with columns `pop`, `region`, `mu`, `sd`, and optional `lb`, `ub` giving the beta’s support (default 0.2 and 1) |
| `h_spec` | Estimation structure for steepness parameters |

The support matters, not just the mean and SD: a beta on $`(0,1)`$ is a
genuinely different function of $`h`$ than one on $`(0.2,1)`$ (it
carries $`\log(h)`$ where the rescaled form carries $`\log(h - 0.2)`$),
and no choice of shape parameters reconciles the two. When bridging an
assessment, match its prior’s support via `lb`/`ub` rather than
approximating with the default.

#### R0 priors

A lognormal prior on $`R_0`$ can be applied per population:

| Argument | Description |
|----|----|
| `use_r0_prior` | 0 = no prior (default), 1 = apply lognormal prior on $`\ln R_0`$ |
| `r0_prior` | Data frame with columns `pop` (population index), `mu` (prior mean on natural scale), and `sd` (prior SD on log scale). Required when `use_r0_prior = 1` |

#### Population straying

When `n_pop > 1`, a fraction of recruits produced by population $`p`$
can “stray” and recruit into regions associated with other populations.

| Argument | Description |
|----|----|
| `stray_rate_spec` | `"fix"` (default), `"est_all"`, etc. |
| `stray_rate_blocks` | Time blocks: `"none_Pop_<p>"` or `"Block_<b>_Year_<s>-<e>_Pop_<p>"` |
| `use_stray_rate_prior` | 0/1 toggle for Beta priors on stray rates |
| `stray_rate_prior` | Data frame with columns `pop`, `block`, `mu` (in (0,1)), `sd` |

#### Spawning movement (single-season, multi-population)

When `n_pop > 1` and `n_seas == 1`, individuals cannot physically move
to their natal region within the seasonal cycle, so a separate
spawning-movement matrix (`sgl_seas_spawning_movement`) routes SSB back
to natal regions for the SRR calculation. If not supplied, SPoRC
defaults to 100% natal homing.

------------------------------------------------------------------------

## Selectivity

Selectivity configuration is shared across fishery fleets
([`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)),
survey fleets
([`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md)),
and retention curves (`ret_sel_model` within
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)).
All three accept the same functional forms and time-varying structures.

#### Functional forms

Specified as character strings following the pattern
`"<form>_Fleet_<f>"` (constant across all time blocks) or
`"<form>_Fleet_<f>_Block_<b>"` (block-specific). The available forms
are:

| String | Functional form | Free parameters |
|----|----|----|
| `"logist1"` | Ascending logistic: $`1 / (1 + \exp(-k(\text{bin} - a_{50})))`$ | 2 ($`\ln a_{50}`$, $`\ln k`$) |
| `"logist2"` | Ascending logistic using $`a_{50}`$ and $`a_{95}`$: $`1 / (1 + 19^{(a_{50} - \text{bin})/a_{95}})`$ | 2 ($`\ln a_{50}`$, $`\ln a_{95}`$) |
| `"gamma"` | Dome-shaped gamma: $`(b/a_\text{max})^{a_\text{max}/p} \exp((a_\text{max} - b)/p)`$ | 2 ($`\ln a_\text{max}`$, $`\ln \delta`$) |
| `"exponential"` | Descending power: $`1/\text{bin}^\beta`$ | 1 ($`\ln \beta`$) |
| `"dbnrml"` | Double-normal with ascending and descending widths, plateau, and endpoint control | 6 |
| `"nonpar"` | Non-parametric: one logit-scale parameter per bin, transformed via $`\text{logit}^{-1}`$ | $`n_\text{bins}`$ |
| `"nonparlog"` | Non-parametric on the log scale, standardized so each year’s selectivity averages to 1 across bins | $`n_\text{bins}`$ |
| `"asymplogist1"` | Logistic with asymptote $`\alpha \in (0,1)`$: $`\alpha / (1 + \exp(-k(\text{bin} - a_{50})))`$ | 3 ($`\text{logit}(\alpha)`$, $`\ln a_{50}`$, $`\ln k`$) |
| `"asymplogist2"` | Logistic with asymptote, $`a_{50}/a_{95}`$ parameterisation | 3 ($`\text{logit}(\alpha)`$, $`\ln a_{50}`$, $`\ln a_{95}`$) |
| `"bicubic"` | Bicubic natural-cubic-spline surface over a bin-node $`\times`$ year-node grid | $`n_\text{bin\_nodes} \times n_\text{yr\_nodes}`$ (see below) |

**Choosing between the two non-parametric forms.** `"nonpar"` bounds
every raw value below one (logistic transform) and mean-standardizes
jointly over years and bins; `"nonparlog"` leaves the scale free and
centers within each year. Under `"nonparlog"` only the within-year
*differences* among parameters are identified, since the level is
absorbed by catchability or fishing mortality, so pair it with the
selectivity parameter centering penalty (below) or fix a bin group via
`*_sel_nonpar_est_bins` to pin the scale. Prefer `"nonpar"` when you
want selectivity interpretable as a proportion without further
constraints.

#### Bicubic spline selectivity (`"bicubic"`)

Unlike the other functional forms above, `"bicubic"` selectivity is
specified with its own extended syntax:

    "bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f>[_Block_<b>][_SelStyr_<year>][_NSelBins_<n>]"

A smooth 2-dimensional selectivity-at-bin-and-year surface is built from
a small grid of $`n_\text{bin\_nodes} \times n_\text{yr\_nodes}`$ freely
estimated log-scale node values (`fish_fixed_sel_pars` /
`srv_fixed_sel_pars`, flattened column-major into a
`[yr_node × bin_node]` matrix). Two natural-cubic-spline weight matrices
are precomputed once at setup:

- $`\mathbf{W}^{\text{bin}}`$
  ($`n_\text{bins} \times n_\text{bin\_nodes}`$), mapping bin-node
  values onto every bin,
- $`\mathbf{W}^{\text{yr}}`$
  ($`n_\text{yrs} \times n_\text{yr\_nodes}`$), mapping year-node values
  onto every year,

and combined via a two-pass tensor product (bin-node values
spline-interpolated across bins for every year-node, then those curves
spline-interpolated across years) to give the full log-selectivity
surface, which is then exponentiated. Setting `n_yr_nodes = 1` collapses
the surface to a time-invariant bin-only spline; combining
`n_yr_nodes = 1` with `fish_sel_blocks`/`srv_sel_blocks` re-fits an
independent bin-only spline within each block.

Two optional suffixes restrict the fitted region of the surface,
edge-holding (flat-lining) outside it:

- `_SelStyr_<year>`: a calendar year within the block. Only years from
  `SelStyr` through the block’s end are actually spline-fit across the
  year dimension; years within the block before `SelStyr` are held
  constant at the `SelStyr` year’s fitted curve (“previous years are
  filled”).
- `_NSelBins_<n>`: restricts the spline fit to the first `n` bins (ages
  or lengths, per `fish_selex_type`/`srv_selex_type`). Bins beyond `n`
  are held constant at the last fitted bin’s value (a plateau), rather
  than continuing the spline extrapolation.

Both suffixes only change which years/bins are considered part of the
*fitted* surface; they do not change the total number of estimated node
parameters.

#### Temporal variation

SPoRC provides two mutually exclusive mechanisms for time-varying
selectivity within a fleet. You can use **discrete blocks** or
**continuous deviations**, but not both on the same fleet.

**Discrete blocks** (`fish_sel_blocks` / `srv_sel_blocks` /
`ret_sel_blocks`): defined as `"Block_<b>_Year_<s>-<e>_Fleet_<f>"` or
`"Block_<b>_Year_<s>-terminal_Fleet_<f>"`. Blocks must be
non-overlapping and collectively span all model years. Each block gets
its own set of fixed-effect selectivity parameters (and can even use a
different functional form).

**Continuous deviations** (`cont_tv_fish_sel` / `cont_tv_srv_sel` /
`cont_tv_ret_sel`): specified as `"<type>_Fleet_<f>"`. For parametric
forms (`"logist1"`, `"gamma"`, etc.), IID and random-walk deviations act
multiplicatively on the transformed base parameters. For semi-parametric
forms (`"3dmarg"`, `"3dcond"`, `"2dar1"`), deviations act
multiplicatively on the selectivity curve at the bin level.

| String | Description |
|----|----|
| `"none"` | Time-invariant (default) |
| `"iid"` | Independent annual deviations on selectivity parameters |
| `"rw"` | Random walk on selectivity parameters |
| `"3dmarg"` | 3D Gaussian Markov random field, marginal variance parameterisation |
| `"3dcond"` | 3D GMRF, conditional variance parameterisation |
| `"2dar1"` | Separable 2D AR(1) over bin × year |

Ancillary controls for continuous time-variation include
`fishsel_pe_pars_spec` (hyperparameter estimation), `fish_sel_devs_spec`
(deviation estimation structure), `fishsel_devs_shared_bins` (bin
grouping for shared deviations), and `corr_opt_semipar` (which
correlation components to suppress in semi-parametric forms).
Bin-sharing specs (`"est_shared_b"` variants) require the deviations to
be indexed by bin: any GMRF/AR1 form, or `"iid"`/`"rw"` on a
non-parametric fleet (`"nonpar"`/`"nonparlog"`, where each deviation
already belongs to one bin); `"iid"`/`"rw"` on a parametric form indexes
deviations by *parameter* and is rejected.

Three further controls tune the process-error penalty itself:

| Argument | Description |
|----|----|
| `fishsel_pe_wt` / `retsel_pe_wt` / `srvsel_pe_wt` | Per-fleet multiplier on the selectivity process-error likelihood (default 1). `0` removes the distributional penalty entirely while the deviations remain estimated; use it when deviations should float subject only to explicit smoothness penalties, as several existing assessments do. Anything other than 0 or 1 makes an estimated PE sigma reinterpretable, so prefer 0/1 unless deliberately down-weighting |
| `fishsel_rw_init_sigma` / `retsel_rw_init_sigma` / `srvsel_rw_init_sigma` | Standard deviation on the first year of an `"rw"` deviation series. Default 5 (first year effectively free). `NA` starts the walk at zero under the walk’s own estimated sigma, making the first year as smooth as every later step; appropriate when the base parametric curve already describes the first year well |
| `*_sel_bin_dev_bins`, `cont_tv_*sel_bin_devs` | Bin-override deviations; see below |

#### Bin-override selectivity deviations

Individual bins can be cut loose from the functional form entirely: bins
named in `fish_sel_bin_dev_bins` / `ret_sel_bin_dev_bins` /
`srv_sel_bin_dev_bins` (a list with one element per fleet, `NULL` for
fleets with no overrides) take a freely estimated annual value
$`\exp(\epsilon_{y,b})`$ in place of whatever the form produced, applied
*after* every other transformation including standardization. The rest
of the curve keeps its parametric shape.

| Argument | Description |
|----|----|
| `fish_sel_bin_dev_bins` etc. | Which bins each fleet overrides (e.g., `list(1, NULL)` frees bin 1 of fleet 1 only) |
| `cont_tv_fishsel_bin_devs` etc. | Process error on the override deviations, per fleet: `"none"`, `"iid"`, or `"rw"` (with its own estimated sigma per bin and `*sel_bin_devs_rw_init_sigma` for the first year) |

**When to use.** The canonical case is a gear whose curve is well
described by a parametric form except for one bin governed by
*availability* rather than the gear (e.g., age-1 availability to a trawl
survey varying with year-class strength). Overriding that bin gives it
free annual variation without abandoning the parametric form (or paying
for a full semi-parametric surface) elsewhere. Give the override `"rw"`
process error unless the bin genuinely jumps independently between
years; with `"none"`, each year’s value is informed only by that year’s
compositions and can be poorly determined in sparse years.

#### Selectivity parameter centering penalty

| Argument | Description |
|----|----|
| `Use_fish_selex_penalty` / `Use_ret_selex_penalty` / `Use_srv_selex_penalty` | 0 (default) / 1 toggle |
| `fish_selex_penalty` etc. | Data frame with columns `region`, `fleet`, `block`, `sex`, `par` (a single index or list-column of integer vectors naming a set), `wt`. Each row penalizes $`w\,[\log(\overline{\exp(\theta)})]^2`$ over the named set |

This pushes the *average selectivity* of a parameter set toward one,
pinning the scalar of a `"nonparlog"` curve that catchability or fishing
mortality would otherwise absorb. It is a softer alternative to fixing a
bin outright, and the standard companion to `"nonparlog"` selectivity
(weights of 10 to 100 are typical starting points; the penalty only
needs to break a ridge, not dominate the fit). Because the expression
averages on the natural scale it is only meaningful for log-scale
parameter sets, so do not apply it to `"nonpar"` (logit-scale) or to
logit-scale asymptote parameters.

#### Parameter sharing and fixing

The `fish_fixed_sel_pars_spec` / `srv_fixed_sel_pars_spec` argument
controls how base selectivity parameters are estimated:

| String             | Meaning                                 |
|--------------------|-----------------------------------------|
| `"est_all"`        | Fully region-, sex-, and fleet-specific |
| `"est_shared_r"`   | Shared across regions                   |
| `"est_shared_s"`   | Shared across sexes                     |
| `"est_shared_r_s"` | Shared across both regions and sexes    |
| `"fix"`            | Fixed at starting values                |

#### Normalisation and length-based selectivity

Selectivity can be normalised relative to a specific bin, to the
maximum, or to the mean across a bin range. When `fit_lengths = 1`,
selectivity operates on length bins and is mapped to age-space via a
user-supplied size-age transition matrix (`SizeAgeTrans`).

#### Selectivity priors

Lognormal priors on selectivity parameters are toggled via
`Use_fish_selex_prior` / `Use_srv_selex_prior`, with hyperparameters
supplied in a data frame (`fish_selex_prior` / `srv_selex_prior`)
containing `region`, `fleet`, `block`, `sex`, `par`, `mu`, and `sd`.

#### Selectivity smoothness penalty weights

All selectivity smoothness/regularisation penalty weights are configured
in a single place,
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md),
via the `fish_sel_pen_wts`, `ret_sel_pen_wts`, and `srv_sel_pen_wts`
arguments (one list per selectivity surface). Each is a named
list/vector; any name not supplied defaults to `0` (off). All six terms
are evaluated directly on a fleet’s *realized*
selectivity-at-bin-at-year surface, so they apply to any selectivity
functional form and any fleet, regardless of block/deviation structure.

| Name | Description |
|----|----|
| `smooth_dome` | Hinge penalty discouraging decreases across adjacent bins (dome-shape control) |
| `smooth_bin_curve` | Second-difference penalty across bins, normalised by the number of fitted bins |
| `smooth_bin_diff` | Unconditional first-difference penalty across bins (both increases *and* decreases contribute, unlike `smooth_dome`), normalised by the number of fitted bins |
| `smooth_yr_diff` | First-difference penalty across years, normalised by the number of fitted years |
| `smooth_yr_curve` | Second-difference penalty across years, normalised by the number of fitted years |
| `smooth_mean_center` | Penalises the per-year mean of log-selectivity away from zero; resolves the scale indeterminacy of the bicubic surface (a uniform per-year shift in log-selectivity otherwise trades off exactly against that year’s fishing mortality) |

**Per-fleet specifications.** A single named specification is shared by
every fleet; alternatively, pass an *unnamed list* with one named
specification per fleet (use
[`list()`](https://rdrr.io/r/base/list.html) or `NULL` for fleets with
no penalties) so, e.g., two surveys can carry different smoothing.

**Per-year weights, bin ranges, and normalization.** Each weight may
also be a vector with one value per model year (`0` skips that year), so
a penalty can act only in years where selectivity changes, or with
year-specific strength. A specification may additionally carry:

| Name | Description |
|----|----|
| `bin_range` | Length-two vector `c(first, last)` restricting the bins the penalties act over; or a named list giving each term its own range. Confine a shape penalty (dome, curvature) to the older ages where the curve should flatten without constraining the ascending limb |
| `normalize` | `TRUE` (default) divides bin-wise weights by the number of penalized bins and year-wise weights by the number of years; settable per term. Turn off when weights are calibrated as explicit variances |
| `yr_diff_ref` | Reference log-selectivity vector for `smooth_yr_diff`’s first penalized year, which otherwise has no predecessor and goes unpenalized. Anchors an otherwise free series to a known selectivity before the data begin |

A useful identity: with `normalize = FALSE`, a per-year `smooth_yr_diff`
weight of $`1/(2\sigma_y^2)`$ is exactly the negative log-kernel of a
random walk on the realized log-selectivity with year-specific standard
deviation $`\sigma_y`$. This is how selectivity random walks with
tabulated per-year sigmas are reproduced as penalties on the curve
rather than as deviation parameters.

**Recommendations.** Prefer the distributional process-error forms
(`cont_tv_*`) when you want an estimated, interpretable sigma; prefer
these penalties when bridging assessments whose selectivity smoothing is
a tuned penalty, or when regularizing a `"bicubic"`/`"nonparlog"`
surface. `smooth_mean_center` (or the centering penalty above) should
accompany any form with a free scale. Weights are on the scale of
squared log-selectivity differences; start small (1 to 10), inspect the
realized surfaces, and remember that normalized and unnormalized weights
differ by a factor of $`n_\text{bins}`$ or $`n_\text{yrs}`$.

------------------------------------------------------------------------

## Catchability

Configured alongside selectivity in
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)
/
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md).

#### Estimation structure (`fish_q_spec` / `srv_q_spec`)

| String           | Description                                        |
|------------------|----------------------------------------------------|
| `"est_all"`      | Free parameter for each region × fleet combination |
| `"est_shared_r"` | One parameter shared across regions (per fleet)    |
| `"fix"`          | Held at user-supplied value                        |

#### Time blocks (`fish_q_blocks` / `srv_q_blocks`)

Same syntax as selectivity blocks: `"none_Fleet_<f>"` for a single
block, or `"Block_<b>_Year_<s>-<e>_Fleet_<f>"` for structured change
points.

#### Analytic catchability (`srv_q_type` / `fish_q_type`)

Per fleet, catchability can be concentrated out of the likelihood rather
than estimated. `srv_q_type` is set in
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md)
and `fish_q_type` in
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md);
the two take the same strings and behave identically, applied to
whichever index that fleet fits:

| String | Description |
|----|----|
| `"est"` | Estimate `ln_srv_q` / `ln_fish_q` (default) |
| `"arith"` | Solve analytically as the ratio of mean observed to mean predicted index |
| `"geo"` | Solve analytically on the log scale, $`\hat{q} = \exp(\overline{\log \text{obs} - \log \text{pred}})`$ |

Both analytic forms solve one $`q`$ per region and fleet from only the
years with observations, ignore block structure, and automatically fix
that fleet’s catchability parameter regardless of `srv_q_spec` /
`fish_q_spec`; they are incompatible with catchability covariates and
priors on that fleet.

**Recommendations.** Use `"geo"` with a lognormal index likelihood: it
is the exact maximum-likelihood solution under a shared SE, removes one
ridge-prone parameter per fleet, and typically speeds and stabilizes
optimization. `"arith"` exists to match assessments that use the
arithmetic ratio. Stay with `"est"` whenever you need a $`q`$ prior,
covariates, time blocks, or when the index’s absolute scale is genuinely
informative (e.g., a swept-area survey with a strong prior near 1, where
concentrating $`q`$ out would discard that information).

#### Priors

| Argument | Description |
|----|----|
| `Use_fish_q_prior` / `Use_srv_q_prior` | 0/1 toggle |
| `fish_q_prior` / `srv_q_prior` | Data frame with `region`, `fleet`, `block`, `mu` (natural scale), `sd` (log scale). Penalty: $`\text{Normal}(\ln(\mu), \sigma)`$ |

------------------------------------------------------------------------

## Natural Mortality

Controlled via
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

#### Estimation vs. fixing (`M_spec`)

| String | Description |
|----|----|
| `"est_ln_M"` | Estimate $`\ln M`$ across the defined block structure |
| `"fix"` | Fix at values supplied via `Fixed_natmort` (array `[p × r × y × a × s]`) |

#### Block structure

Natural mortality parameters can be shared or made specific along any
combination of five axes. Each argument takes either `'constant'` (all
levels pooled) or a list of integer vectors defining blocks:

| Argument | Axis | Example |
|----|----|----|
| `M_popblk_spec` | Population | `list(c(1,2), 3)` → pops 1-2 share $`M`$, pop 3 is separate |
| `M_regionblk_spec` | Region | `list(1:3, 4:5)` |
| `M_yearblk_spec` | Year | `list(1:20, 21:40)` |
| `M_ageblk_spec` | Age | `list(1:5, 6:30)` → young vs. old |
| `M_sexblk_spec` | Sex | `list(1, 2)` → sex-specific $`M`$ |

All block specifications are crossed to produce unique $`M`$ parameters.
For instance, two age blocks × two sex blocks = four estimated $`\ln M`$
values (assuming everything else is `'constant'`).

#### Priors

Lognormal priors on $`M`$ are activated via `Use_M_prior = 1`, supplying
a data frame with columns `popblk`, `regionblk`, `yearblk`, `ageblk`,
`sexblk`, `mu`, and `sd`.

------------------------------------------------------------------------

## Movement

Configured via
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md)
(EM) or `Setup_Sim_Movement()` (OM).

SPoRC implements two movement parameterisations, plus a fixed-matrix
escape hatch.

#### Movement type (`move_type`)

| Value | Description |
|----|----|
| 0 | **Unstructured Markov.** Transition probabilities estimated via multinomial logit (softmax), with region 1 as the implicit reference category. Supports blocking along population, age, year, season, and sex dimensions |
| 1 | **CTMC.** A continuous-time Markov chain builds an instantaneous rate matrix $`Q = D + Z`$ from diffusion ($`D`$, isotropic dispersal scaled by region area) and taxis ($`Z`$, directional preference). Transition probabilities are obtained by matrix exponentiation: $`\mathbf{P} = \exp(Q)`$ |

#### Fixed movement (`use_fixed_movement = 1`)

Bypasses estimation entirely. Supply `Fixed_Movement` as an array
`[p × r_\text{from} × r_\text{to} × y × \tau × a × s]`.

#### Unstructured Markov block structure (`move_type = 0`)

Each dimension can be pooled (`'constant'`) or blocked:

| Argument | Description |
|----|----|
| `Movement_popblk_spec` | Population blocks |
| `Movement_ageblk_spec` | Age blocks (e.g., `list(1:5, 6:30)` for age-dependent movement) |
| `Movement_yearblk_spec` | Year blocks |
| `Movement_seasblk_spec` | Season blocks |
| `Movement_sexblk_spec` | Sex blocks |

#### CTMC configuration (`move_type = 1`)

| Argument | Description |
|----|----|
| `ctmc_move_dat` | Data frame of covariates with required columns `pop`, `regions`, `years`, `seas`, `ages`, `sexes`, plus any variables referenced in the formulas. For projection years, covariate lookups are capped at the last historical year unless extended rows are supplied |
| `diffusion_formula` | R formula for diffusion (e.g., `~ 1`, `~ bs(years, df = 4)`) |
| `preference_formula` | R formula for taxis/preference |
| `adjacency_mat` | Square `[n_regions × n_regions]` binary connectivity matrix |
| `area_r` | Numeric vector of region areas (scales diffusion rates) |
| `ctmc_diffusion_bounds` | 0/1: if `1`, shifts diffusion columns to guarantee all off-diagonal generator-matrix entries are non-negative |
| `ctmc_scale_by_seasdur` | 0/1 (default `1`): if `1`, treats the generator as an annual rate and exponentiates `Q * seasdur[s]` each season, so movement and mortality share time units. `0` results in `Q * 1` each season. |

#### Movement and mortality sequencing (`move_timing`)

Controls how movement and mortality are ordered within a season. The
three options coincide exactly when total mortality is constant across
regions, when mortality is zero, and when movement is absent, so
single-region models are unaffected by the setting.

| Value | Sequencing |
|----|----|
| `0` (default) | Movement then mortality. Historical `SPoRC` behaviour |
| `1` | Mortality then movement |
| `2` | Continuous: movement and mortality act simultaneously, via `expm(Q * seasdur - diag(Z))` |

`move_timing = 2` requires an estimated CTMC generator (`move_type = 1`
with `use_fixed_movement = 0`); a discrete movement matrix has no
guaranteed real generator, so it is rejected rather than approximated.

The setting is not confined to the projection step. Anything observed
partway through a season has to be evaluated consistently with it:

| Quantity | `0` | `1` | `2` |
|----|----|----|----|
| Spawning biomass | Post-movement location | Pre-movement location | Partially redistributed, `expm(A * t_spawn)` |
| Catch, discards | Baranov, post-movement | Baranov, pre-movement | Season-integrated abundance (spatial Baranov) |
| Fishery index | Post-movement `N` | Pre-movement `N` | Season-integrated abundance |
| Survey index | `N * exp(-t_srv * Z)` | `N * exp(-t_srv * Z)` | Snapshot, `expm(A * t_srv) N` |
| Tag recaptures | Baranov on tag cohort | Baranov on tag cohort | Season-integrated tag abundance |

where `A = t(Q) * seasdur - diag(Z)`. Equilibrium initialisation
(including the matrix plus-group solution), per-recruit reference
points, and projections all use the same seasonal operator, so they
inherit the setting as well. Note the survey/fishery contrast: a survey
is a snapshot at an instant within the season and so uses a partial
propagation, whereas catch and the fishery index accumulate over the
season and so use the integral of that propagation. Full equations are
in the *Description of Model Equations* vignette.

#### Continuous movement deviations (`cont_vary_movement`)

Origin-destination deviations applied multiplicatively to off-diagonal
rates (CTMC) or additively on the logit scale (unstructured).

| String               | Deviation structure                    |
|----------------------|----------------------------------------|
| `"none"`             | No deviations (default)                |
| `"iid_y"`            | Year only                              |
| `"iid_a"`            | Age only                               |
| `"iid_y_a"`          | Year × age                             |
| `"iid_y_a_s"`        | Year × age × sex                       |
| `"iid_y_seas_a_s"`   | Year × season × age × sex              |
| `"iid_p_y"`          | Population × year                      |
| `"iid_p_a"`          | Population × age                       |
| `"iid_p_y_a"`        | Population × year × age                |
| `"iid_p_y_a_s"`      | Population × year × age × sex          |
| `"iid_p_y_seas_a_s"` | Population × year × season × age × sex |

#### Additional movement controls

| Argument | Description |
|----|----|
| `do_recruits_move` | 0 = age-1 fish do not move; 1 = recruits follow the movement matrix |
| `Use_Movement_Prior` | 0/1 toggle for Dirichlet movement priors (unstructured) |
| `Movement_prior` | Data frame with `pop`, `region_from`, `year`, `seas`, `age`, `sex`, `alpha` (Dirichlet concentration vector of length `n_regions`) |
| `Movement_cont_pe_pars_spec` | Estimation structure for process-error hyperparameters: `"none"`, `"fix"`, `"est_all"`, `"est_shared"` |

------------------------------------------------------------------------

## Catch, Fishing Mortality, and Discards

Configured via
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

#### Catch conditioning

| Argument | Description |
|----|----|
| `ObsCatch` | Observed aggregate catch `[r × y × τ × f]` |
| `ObsCatch_pop` | Population-specific catch `[p × r × y × τ × f]` |
| `catch_units` | Per-fleet: 0 = abundance, 1 = biomass |
| `Use_F_pen` | 0/1 toggle for the fishing mortality deviation penalty (`Fmort_nLL`); see [Which deviations are penalized](#which-deviations-are-penalized) |
| `sigmaC_spec` / `sigmaC_pop_spec` | Catch observation-error estimation |

#### Fishing mortality process error

`ln_F_devs` (annual log-scale deviations about `ln_F_mean`) can follow
one of three process-error structures, set via `Fdev_model`:

| String  | Description                             |
|---------|-----------------------------------------|
| `"iid"` | Independent annual deviations (default) |
| `"rw"`  | Random walk                             |
| `"ar1"` | First-order autoregressive              |

**`Fdev_pen_center`** (`"fixed"`, the default, or `"own_mean"`) sets
where the iid deviation penalty is centred; see the deviation penalty
centering discussion in the Recruitment section for the shared
rationale. Under the mean-plus-deviations parameterization the level of
$`F`$ is already carried by `ln_F_mean`, so `"own_mean"` penalizes only
the deviations’ spread and avoids constraining the level twice, matching
assessments whose deviation vectors sum to zero. It leaves `ln_F_mean`
and the deviations’ level mutually unidentified unless one is fixed,
which `ln_F_mean_spec = "fix"` does; setup warns when it is combined
with an estimated mean in configurations where nothing else reads the
mean. Note also that the default iid penalty with `"fixed"` centering
acts as a hidden “biomass tracks catch” prior when catch trends
strongly; `Fdev_model = "rw"` is usually the better remedy for that than
re-centering.

**`ln_F_mean_spec`** (`"est"`, the default, or `"fix"`) chooses the
fishing mortality parameterization. `"est"` is the mean-plus-deviations
form. `"fix"` maps `ln_F_mean` off at its starting value (zero unless
supplied), so the deviations are free annual log-F outright:
$`F = \exp(\epsilon)`$. It must be paired with
`Fdev_pen_center = "own_mean"`, `Fdev_model = "rw"`, or `Use_F_pen = 0`
— an `"iid"` or `"ar1"` penalty centred on the fixed zero mean would
shrink the deviations toward $`F = 1`$, and setup rejects that
combination.

`sigmaF_spec` controls sharing/fixing of the process-error standard
deviation (`ln_sigmaF`) across region × season × fleet, using the same
`"est_all"` / `"est_shared_<dims>"` / `"fix"` convention described in
[Parameter sharing conventions](#parameter-sharing-conventions).
`Fdev_rho_spec` analogously controls the AR1 correlation parameter
(`Fdev_rho`) and is only active when `Fdev_model = "ar1"`; for `"iid"`
or `"rw"`, `Fdev_rho` is unused and mapped to `NA` regardless of what is
supplied.

Catch-active years (`UseCatch == 1` or any `UseCatch_pop == 1`) do
**not** need to be contiguous under `"rw"` or `"ar1"`, a fishery may
close for several years and reopen later. The transition between two
active years separated by a gap of $`d`$ closed years is taken directly
over the elapsed gap (the same marginal transition obtained by
estimating deviations for the closed years and integrating them out,
without actually estimating them). See \[Fishing Mortality Deviations\]
in
[`vignette("c_model_equations")`](https://chengmatt.github.io/SPoRC/dev/articles/c_model_equations.md)
for the exact equations, and [Which deviations are
penalized](#which-deviations-are-penalized) for how the active sequence
is determined.

#### Discard and retention framework

SPoRC decomposes total fishing mortality at age into retained and
dead-discard components:

``` math
F_{p,r,\tau,a,s} = \sum_f F_{r,\tau,f} \left[\underbrace{\text{sel}_{a,f} \cdot \text{ret}_{a,f}}_{\text{retained}} + \underbrace{\text{sel}_{a,f} \cdot (1 - \text{ret}_{a,f}) \cdot \text{dmr}_f}_{\text{dead discards}}\right]
```

| Argument | Description |
|----|----|
| `ret_sel_model` | Retention selectivity functional form (same options as `fish_sel_model`) |
| `ret_sel_blocks` | Retention selectivity time blocks |
| `cont_tv_ret_sel` | Continuous time-varying retention selectivity |
| `dmr_mean_spec` | Estimation structure for dead discard mortality rate means |
| `dmr_dev_spec` | Estimation structure for DMR deviations: `"fix"` (default) or `"est_all"` |
| `Use_dmr_pen` | 0/1 toggle for the DMR deviation penalty (`dmr_nLL`); must be `1` when `dmr_dev_spec = "est_all"` |
| `discard_units` | Per-fleet: 0 = abundance, 1 = biomass, 2 = abundance fraction, 3 = biomass fraction |
| `ObsDiscard` / `ObsDiscard_pop` | Observed aggregate and population-specific discards |
| `sigmaD_spec` / `sigmaD_pop_spec` | Discard observation-error estimation |

#### Estimating discard mortality rates

`dmr_dev_spec = "est_all"` estimates an annual deviation in every
**fished** region × year × season × fleet cell, not only in cells
carrying a discard observation. This is deliberate. From the
decomposition above, `dmr` scales the dead-discard component of total
mortality, so it enters $`Z`$ and therefore propagates into predicted
retained catch, indices, and compositions:

``` math
C_{a} = \frac{F^{\text{ret}}_{a}}{Z_{a}} N_{a} \left(1 - e^{-Z_{a}}\right)
```

A cell is therefore informative about `dmr` whenever it is fished and
retention is less than one, with or without discard data. Discard
observations are not the dividing line, `dmr` in fact cancels out of the
predicted discard itself, which divides the dead discards back through
by it, leaving only the same $`Z`$ dependence.

The cells excluded are true closures, matching the condition under which
the objective pins `dmr` to zero: no aggregate or population-specific
catch is fit, and the aggregate catch observation is a recorded zero
rather than missing (`NA`). A missing observation is treated as a fished
year, on the assumption that fishing continued and we simply lack a
value to fit.

Note that a fleet with full retention (`ret_sel = 1`, the package
default) contributes no dead discards at all, so `dmr` drops out of the
objective entirely and its deviations have a flat gradient.
`dmr_dev_spec` and `dmr_mean_spec` are only meaningful once retention is
actually modelled, see
[`vignette("s_discard_retention")`](https://chengmatt.github.io/SPoRC/dev/articles/s_discard_retention.md).

#### Which deviations are penalized

The `ln_F_devs`, `logit_dmr_devs` and `ln_RecDevs` penalties are
evaluated on exactly the deviations that are estimated, read off the map
mirrors `map_ln_F_devs`, `map_logit_dmr_devs` and `map_ln_RecDevs` in
the data list rather than recomputed from the catch indicators or the
deviation index ranges.
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
refreshes those mirrors from the map immediately before building the
model, so a deviation mapped off by hand after setup is neither
estimated nor penalized:

``` r

# drop dmr deviations for a stretch of years, penalty included
mp <- array(input_list$map$logit_dmr_devs, dim = dim(input_list$par$logit_dmr_devs))
mp[, 10:20, , ] <- NA
input_list$map$logit_dmr_devs <- factor(mp)
```

Two consequences worth knowing. Mapping a deviation off pins it at
whatever sits in `$par`, which is `0` by default, so `dmr` falls back on
`plogis(logit_dmr_mean)`, but a starting value supplied through `...` is
kept, not reset to zero. And under `Fdev_model = "rw"` or `"ar1"`, a
pinned deviation is dropped from the active sequence, widening the gap
$`d`$ between the deviations either side of it rather than being treated
as an actual year.

For `ln_RecDevs` this matters most when `ln_sigmaR` is estimated. A
deviation fixed at zero still sits at the penalty’s mean, so were it
penalized it would add to the sum of squares’ denominator without adding
any spread, and $`\sigma_R`$ would be pulled low by roughly
$`\sqrt{n_{est} / n_{total}}`$. Holding the first 80 of 255 quarterly
deviations, for instance, biases $`\sigma_R`$ down by about 17%.

With `ln_sigmaR` fixed the excluded terms are constants: the objective
shifts by a fixed amount and its gradient is unchanged everywhere, so
the likelihood surface and the location of any optimum are untouched,
and `Rec_nLL` becomes comparable across runs that hold different numbers
of deviations. Note that this is not the same as the *fitted values*
being unchanged. `nlminb` tests convergence on relative changes in the
objective’s value, so on a model that stops short of convergence (a flat
ridge, or a `"singular convergence"` or `"false convergence"` message),
shifting the objective by a constant can move where the optimizer halts.
If excluding these terms visibly changes your estimates, the
identifiability of the model is the thing to look at, not the penalty.

------------------------------------------------------------------------

## Biological Inputs

Supplied via
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md)
(EM) or
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md)
(OM).

| Input | Dimensions | Description |
|----|----|----|
| `WAA` | `[p × r × y × τ × a × s]` | Spawning weight-at-age (used for SSB) |
| `WAA_fish` | `[p × r × y × τ × a × s × f]` | Fishery-specific weight-at-age. Defaults to `WAA` if `NULL` |
| `WAA_srv` | `[p × r × y × τ × a × s × f]` | Survey-specific weight-at-age. Defaults to `WAA` if `NULL` |
| `MatAA` | `[p × r × y × τ × a × s]` | Maturity-at-age (proportions in \[0, 1\]) |
| `AgeingError` | `[a_\text{model} × a_\text{obs}]` or `[y × a_\text{model} × a_\text{obs}]` | Row-stochastic matrix mapping true ages to observed bins. Defaults to identity. Supply explicitly when observed ages are a subset of modelled ages |
| `SizeAgeTrans` | `[p × r × y × τ × l × a × s]` | Size-age transition matrix. Required when `fit_lengths = 1` |
| `fit_lengths` | 0/1 | Toggle for fitting length compositions |
| `addtocomp` | scalar | Small constant guarding $`\log(0)`$ in composition likelihoods (default `1e-3`) |
| `comp_const_obs` | 0/1 | Whether `addtocomp` is also added to the observed proportions used as multinomial weights. `1` (default) is the long-standing SPoRC behaviour; `0` adds it only inside the logarithms, a convention several existing assessments use. The difference slightly reweights every bin, so set `0` when bridging such assessments and otherwise leave the default |

------------------------------------------------------------------------

## Observation Model: Indices and Compositions

Indices and compositions are configured together for fishery fleets
([`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md))
and survey fleets
([`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)).

#### Indices of abundance

| Argument | Description |
|----|----|
| `ObsFishIdx` / `ObsSrvIdx` | Observed index values |
| `ObsFishIdx_SE` / `ObsSrvIdx_SE` | Log-scale standard errors |
| `fish_idx_type` / `srv_idx_type` | 0 = abundance, 1 = biomass |
| `ObsFishIdx_pop` / `ObsSrvIdx_pop` | Population-specific indices (separate likelihood contribution) |
| `ObsFishIdx_pop_SE` / `ObsSrvIdx_pop_SE` | Population-specific index SEs |

#### Index error structure (`FishIdx_LikeType` / `SrvIdx_LikeType`)

Per fleet. A fleet’s population-specific index stream follows the same
choice for `"lognormal"` and `"normal"`, but stays lognormal under
`"mvn"`, whose covariance describes the region-aggregated series only.
The simulator (`Setup_Sim_Fishing` / `Setup_Sim_Survey`, and
`simulation_self_test` automatically) draws index observations under the
same structures, with an `"mvn"` fleet drawn through a one-factor
decomposition of its covariance:

| String | Description |
|----|----|
| `"lognormal"` | Default. SEs are log-scale standard deviations |
| `"normal"` | Normal on the arithmetic scale; SEs are arithmetic standard deviations |
| `"mvn"` | Multivariate normal on the arithmetic scale with a fixed covariance supplied via `FishIdx_Cov` / `SrvIdx_Cov` (one matrix per `"mvn"` fleet, square with one row per fitted observation, ordered as observations appear scanning the fleet’s use flags in array order). Validated at setup for symmetry and positive definiteness |

**Recommendations.** `"lognormal"` is right for almost all abundance
indices (positive, multiplicative errors). Use `"mvn"` when the index
provider supplies a covariance across years (e.g., a model-based index
from VAST or sdmTMB whose inter-annual correlations are real information
a diagonal likelihood would double-count). Use `"normal"` only for
series that can legitimately go near zero or negative, or when bridging
an assessment that fits on the arithmetic scale. Two caveats: OSA
residuals are only available for `"lognormal"` fleets, and the MVN
density includes the $`-\tfrac{n}{2}\log(2\pi)`$ constant that some
other implementations omit, so absolute likelihood values are not
directly comparable across implementations even when fits match.

#### Index timing and age restriction

| Argument | Description |
|----|----|
| `t_fish` | `[r × τ × f]`, fraction of the season elapsed when the fishery index is observed; numbers decay by $`e^{-t Z}`$ before the index is formed, mirroring `t_srv`. Default `0` (start of season, the historical behaviour). Set `0.5` for a mid-season CPUE snapshot |
| `fish_idx_ages` / `srv_idx_ages` | Which ages contribute to each fleet’s index total: a list (one vector of ages per fleet, `NULL` = all) or an `[n_ages × n_fleets]` 0/1 array. The restriction applies to the index *sum* only; selectivity, catch, and the fleet’s compositions are untouched. Restricting a survey to one age turns it into an index of that age alone (e.g., an age-1 acoustic recruitment index), without needing a knife-edge selectivity that would corrupt the compositions |

#### Composition bin restriction (`FishAgeComps_bins` / `SrvAgeComps_bins`)

Same list-or-array format as `*_idx_ages`, but restricting which
*observed* bins (i.e., after ageing error) a fleet’s age compositions
are fit over. Observed and expected compositions are both subset and
renormalized within the named bins; bins outside are simply left out of
the likelihood rather than forced to be explained. Use for gears that
never resolve part of the age range (e.g., a fishery that never catches
the youngest ages, whose structural zeros would otherwise carry
information); leave alone when the zeros are informative sampling zeros.

#### Composition likelihood families

Age and length compositions each accept one of five likelihood families:

| Value | Likelihood |
|----|----|
| 0 | Multinomial |
| 1 | Dirichlet-multinomial (overdispersion parameter $`\theta`$ estimated per fleet) |
| 2 | Logistic-normal, independent bins |
| 3 | Logistic-normal with AR(1) correlation across bins |
| 4 | Logistic-normal with AR(1) bin correlation and constant cross-sex correlation |

These are set per data stream via `comp_fishage_like`,
`comp_fishlen_like`, `comp_srvage_like`, `comp_srvlen_like`, and their
`_pop` and `_discard` variants.

#### Composition structure types

Each composition data stream has a “type” controlling how data are
aggregated:

| Value | Structure |
|----|----|
| 0 | Aggregated across sexes and regions |
| 1 | Split by sex and region (no implicit sex-ratio information) |
| 2 | Joint across sexes, split by region (preserves sex-ratio information) |
| 999 | No data for this fleet/year |

#### Data streams

SPoRC supports a full matrix of composition data streams, each
independently configurable:

| Category | Aggregate | Population-specific |
|----|----|----|
| Fishery age | `comp_fishage_like` | `comp_fishage_pop_like` |
| Fishery length | `comp_fishlen_like` | `comp_fishlen_pop_like` |
| Fishery discard age | `comp_fishage_discard_like` | `comp_fishage_discard_pop_like` |
| Fishery discard length | `comp_fishlen_discard_like` | `comp_fishlen_discard_pop_like` |
| Survey age | `comp_srvage_like` | `comp_srvage_pop_like` |
| Survey length | `comp_srvlen_like` | `comp_srvlen_pop_like` |

Each stream carries its own ISS arrays, $`\theta`$ parameters (log-scale
overdispersion), and correlation parameters.

------------------------------------------------------------------------

## Tagging

Configured via
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
(EM) or
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md)
(OM).

SPoRC implements a conventional mark-recapture framework (Brownie-type
likelihood) where tagged individuals follow the full seasonal population
dynamics, movement, mortality, and fishing, after release.

| Feature | Description |
|----|----|
| Release structure | Tags released by population, region, season, age, and sex |
| Recapture structure | Fleet-specific recaptures with release-cohort × recapture-year likelihood |
| Tag shedding | Chronic tag-loss rate (estimated or fixed) |
| Initial tag mortality | Immediate post-release mortality parameter |
| Reporting rates | Fleet-specific tag-reporting rates (estimated or fixed) |
| Population/age/sex attribution | Tags can carry population, age, and sex information at release |

------------------------------------------------------------------------

## Reference Points

Computed post-estimation via
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md).

The function accepts a `type` argument for spatial structure and a
`what` argument for the reference-point method. All methods project a
single recruit through the full age, season, and spatial structure to
compute SBPR and YPR.

#### Available methods

| `type` | `what` | Description | Multi-pop? |
|----|----|----|----|
| `"single_region"` | `"SPR"` | $`F_{\text{SPR}_x}`$, no movement | , |
| `"single_region"` | `"BH_MSY"` | Beverton-Holt $`F_\text{MSY}`$, no movement | , |
| `"multi_region"` | `"independent_SPR"` | Per-region $`F_{\text{SPR}_x}`$ ignoring movement | ✓ |
| `"multi_region"` | `"independent_BH_MSY"` | Per-region $`F_\text{MSY}`$ ignoring movement | ✓ |
| `"multi_region"` | `"global_SPR"` | Single $`F_{\text{SPR}_x}`$ applied uniformly across regions, with movement | ✓ |
| `"multi_region"` | `"global_BH_MSY"` | Single $`F_\text{MSY}`$ with movement (single-pop only) | , |
| `"multi_region"` | `"local_BH_MSY"` | Region-specific $`F_\text{MSY}`$ values jointly maximising total yield under movement. Uses Newton-Raphson to solve equilibrium recruitment by origin | ✓ |

#### Controls

| Argument | Description |
|----|----|
| `SPR_x` | Target SPR fraction (e.g., 0.4 for $`B_{40\%}`$) |
| `n_avg_yrs` | Terminal years to average demographic rates (selectivity, $`M`$, WAA, maturity, movement) |
| `calc_rec_st_yr` | First year for computing mean historical recruitment |
| `rec_age` | Recruitment lag used to exclude most-recent years from the mean |
| `is_discard_fleet` | Integer vector (per fleet) flagging discard-only fleets to exclude from landed yield in MSY calculations |
| `local_bh_msy_newton_steps` | Newton iterations for the local $`F_\text{MSY}`$ solve (default 6) |

------------------------------------------------------------------------

## Simulation and Closed-Loop MSE

The operating-model side uses a parallel set of `Setup_Sim_*` functions
that mirror the estimation-model pipeline but populate a `sim_list`
rather than `input_list`. Key simulation-specific functions:

| Function | Purpose |
|----|----|
| [`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md) | Forward-project the OM without feedback (open-loop) |
| [`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/dev/reference/condition_closed_loop_simulations.md) | Closed-loop MSE: periodically re-fits the EM, applies an HCR, and updates $`F`$ |
| [`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_data_to_SPoRC.md) | Converts OM output (with observation error) to `input_list` format for the EM |
| [`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md) | Simulation-estimation test: fits the EM back to OM-generated data |
| [`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/dev/reference/get_closed_loop_reference_points.md) | Computes reference points inside the MSE feedback loop |

Key closed-loop arguments in
[`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/dev/reference/condition_closed_loop_simulations.md):

| Argument | Description |
|----|----|
| `closed_loop_yrs` | Year range over which feedback is active |
| `assessment_period` | Frequency of EM re-fitting (e.g., every 2 years) |
| `use_true_values` | If `TRUE`, the HCR uses OM-truth instead of EM estimates (perfect-information benchmark) |

------------------------------------------------------------------------

## Estimation and Optimisation

[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
constructs the RTMB automatic-differentiation function, optimises via
`nlminb`, and refines with Newton steps.

| Argument | Default | Description |
|----|----|----|
| `data` | , | `input_list$data` |
| `parameters` | , | `input_list$par` |
| `mapping` | , | `input_list$map` |
| `random` | `NULL` | Parameter names to marginalise as random effects (Laplace approximation) |
| `newton_loops` | `3` | Post-convergence Newton steps ($`\Delta\theta = -H^{-1}g`$) to reduce residual gradients |
| `do_optim` | `TRUE` | `FALSE` returns the un-optimised `MakeADFun` object for debugging |
| `nlminb_control` | `list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)` | Passed to [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html) |

------------------------------------------------------------------------

## Diagnostics

| Function | What it does |
|----|----|
| [`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md) | Sequentially peels terminal years and re-fits; reports Mohn’s $`\rho`$ per quantity |
| [`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md) | Refits from perturbed starting values to test convergence stability |
| [`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md) | Profiles the likelihood surface over user-specified parameters |
| [`do_francis_reweighting()`](https://chengmatt.github.io/SPoRC/dev/reference/do_francis_reweighting.md) | Computes Francis TA1.8 weights for composition data |
| [`run_francis()`](https://chengmatt.github.io/SPoRC/dev/reference/run_francis.md) | Iterative Francis reweighting loop (re-fits after each adjustment) |
| [`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md) | One-step-ahead residuals, compositions, conventional tagging, and Catch/Discard/FishIdx/SrvIdx indices (external post-hoc, or internal model-based via `do_internal_comp_osa`/`do_internal_conv_tag_osa`). See [`vignette("u_osa_residuals")`](https://chengmatt.github.io/SPoRC/dev/articles/u_osa_residuals.md) |
| [`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md) | Runs test for serial correlation in residuals |
| [`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md) | Extracts model report quantities across MCMC posterior draws (compatible with `adnuts` / `tmbstan`) |
| [`marg_AIC()`](https://chengmatt.github.io/SPoRC/dev/reference/marg_AIC.md) | Marginal AIC for models with random effects |

------------------------------------------------------------------------

## Plotting

| Function | Output |
|----|----|
| [`plot_all_basic()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_all_basic.md) | Multi-panel diagnostic overview |
| [`get_ts_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_ts_plot.md) | Time series: SSB, recruitment, $`F`$, catch, depletion |
| [`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md) | Index fits (observed vs. predicted) |
| [`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md) | Catch fits |
| [`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md) | Composition fits (bubble plots / proportion plots) |
| [`get_selex_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_plot.md) | Selectivity-at-age or selectivity-at-length |
| [`get_biological_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_biological_plot.md) | WAA, maturity, natural mortality |
| [`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md) | Likelihood component breakdown |
| [`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md) | Retrospective trajectories |
| [`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_data_fitted_plot.md) | Data-availability timeline |
| [`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md) | Residual diagnostics |
| [`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md) | Relative difference plots for retrospective analyses |

All plotting functions accept either a single fitted model or a list of
models for side-by-side comparison.
