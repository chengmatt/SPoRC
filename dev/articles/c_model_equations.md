# Description of Model Equations

The Stochastic Population over Regional Components (`SPoRC`) model is a
generalized integrated population model written in `RTMB` (R bindings
for Template Model Builder; Kristensen et al., 2016) that supports age,
sex, population, seasonal, and spatially-structured dynamics. Population
dynamics operate across an annual time step that is further subdivided
into $`n_\tau`$ seasons of duration $`\Delta\tau`$ (where
$`\sum_\tau \Delta\tau = 1`$). Within each annual-seasonal cycle,
processes occur in the following order:

1.  Recruitment generally occurs in the first season, with additional
    recruits apportioned to subsequent seasons according to seasonal
    proportions. The exception is age-0 recruitment (no lag between
    spawning and recruitment; see Recruitment Processes below), where
    recruits instead enter no earlier than the spawning season,
2.  Markovian movement of individuals (movement only occurs in the
    spatial model) and total mortality both act within the season; at
    the end of the final season, individuals advance in age.

How movement and mortality are sequenced within step 2 is a user choice,
set by `move_timing`: movement then mortality (`0`, the default and the
historical `SPoRC` behavior), mortality then movement (`1`), or the two
acting simultaneously and continuously (`2`). This choice is not
confined to the projection equations. It propagates to every quantity
whose value depends on where fish are partway through a season, spawning
biomass, catch-at-age, fishery and survey indices, tag recaptures,
equilibrium initialization, and per-recruit reference points, and each
of those is given below. The three options coincide exactly when total
mortality is constant across regions, when mortality is zero, and when
movement is absent, so single-region models are unaffected by the
setting.

Tag releases occur simultaneous to recruitment in the release year and
season (i.e., recruits can be tagged), and tag recaptures are computed
each season.

These processes are modeled across five primary partitions: population
($`p`$), region ($`r`$), year ($`y`$), season ($`\tau`$), age ($`a`$ and
$`a_{+}`$, where $`a_{+}`$ is the plus group), and sex ($`s`$). In
single-population, single-region and/or single-sex models, these
equations generally reduce by setting $`p = 1`$, $`r = 1`$ and/or
$`s = 1`$. In general, the same equations are used for both simulation
and estimation.

### Process Equations

#### Population Initialization

In `SPoRC`, three primary methods exist to initialize the equilibrium
population of the model. The first method derives the equilibrium
population using the following process:

``` math
N_{p,r,a,s}^{'} = \mu_p^{\text{RecInit}}\exp\left( - (a - 1) \cdot Z_{p,r,a,s}^{'} \right)\psi_{p,r,y = 1,s}\zeta_{p,r},\quad\text{for }2 \leq a < a_{+}
```

``` math
Z_{p,r,a,s}^{'} = \text{Natmort}_{p,r,y = 1,a,s} + \sum_{f=1}^{n_f} F_{r,\tau,f}^{\text{Init}} \cdot \mathbb{1}_{r,\tau,f}^{\text{Catch}} \cdot \left[\text{Sel}_{p,r,y=1,\tau,a,s,f}^{\text{Fsh}} \cdot \text{Sel}_{p,r,y=1,\tau,a,s,f}^{\text{Ret}} + \text{Sel}_{p,r,y=1,\tau,a,s,f}^{\text{Fsh}} \cdot \left(1 - \text{Sel}_{p,r,y=1,\tau,a,s,f}^{\text{Ret}}\right) \cdot \delta_{r,\tau,f}\right]
```

The initialization fishing mortality $`F_{r,\tau,f}^{\text{Init}}`$ is
derived from a single parameter $`\theta_{r,\tau,f}^{\text{Init}}`$
(`init_F_par`), in one of two forms set by `init_F_form`:

``` math
F_{r,\tau,f}^{\text{Init}} = \begin{cases}
\text{logit}^{-1}\left(\theta_{r,\tau,f}^{\text{Init}}\right) \cdot \exp\left(\mu_{r,\tau,f}^{\text{Fsh}}\right), & \texttt{init\_F\_form = "prop"}\\[4pt]
\exp\left(\theta_{r,\tau,f}^{\text{Init}}\right), & \texttt{init\_F\_form = "abs"}
\end{cases}
```

Under `"prop"` the initialization F is a proportion of the mean fishing
mortality, bounded to $`(0,1)`$ by the inverse-logit, so the initial age
structure moves with $`\mu_{r,\tau,f}^{\text{Fsh}}`$. Under `"abs"` it
is an absolute rate and is independent of
$`\mu_{r,\tau,f}^{\text{Fsh}}`$. Whether $`\theta^{\text{Init}}`$ is
estimated is set separately by `init_F_spec` (`"fix"` or `"est"`), so
all four combinations are available.

The distinction matters whenever the initial condition is fished. Catch
constrains only the product of numbers-at-age and fishing mortality, so
under `"prop"` a single parameter both depletes the initial age
structure and scales the F series: the optimizer can fit the observed
catch equally well with a smaller, harder-fished stock. Use `"abs"` when
the historical fishing mortality that shaped the initial condition is
conceptually distinct from the mean F of the modeled period, for example
when bridging an assessment that carries its own separate historical F
parameter.

where:

- $`N_{p,r,a,s}^{'}`$ are the equilibrium numbers-at-age,
- $`\mu_p^{\text{RecInit}}`$ is a global recruitment parameter used to
  scale the equilibrium age structure during initialization. Users have
  the option to either initialize the population using the same
  recruitment parameter that governs the stock-recruit relationship
  (either virgin or mean recruitment depending on the parameterization),
  or to estimate a separate recruitment scalar exclusively for
  initialization The latter is useful when the historical mean
  recruitment used to initialize the population differs from the virgin
  recruitment implied by the stock-recruit relationship, or when the
  assumption that the population was at unfished equilibrium at the
  start of the time series is not appropriate.
- $`Z_{p,r,a,s}^{'}`$ is the initial instantaneous total mortality rate,
- $`\text{Natmort}_{p,r,y = 1,a,s}`$ is the instantaneous natural
  mortality rate,
- $`\mu_{r,\tau,f}^{\text{Fsh}}`$ is the log-mean fishing mortality rate
  for fleet $`f`$ in region $`r`$ and season $`\tau`$,
- $`F_{r,\tau,f}^{\text{Init}}`$ is the initialization fishing mortality
  for fleet $`f`$, either a proportion of the mean fishing mortality or
  an absolute rate (see below),
- $`\mathbb{1}_{r,\tau,f}^{\text{Catch}}`$ is an indicator variable
  equal to 1 if fleet $`f`$ is active in region $`r`$ and season
  $`\tau`$ in year 1, and 0 otherwise,
- $`\text{Sel}_{p,r,y = 1,\tau,a,s,f}^{\text{Fsh}}`$ is the total
  fishery selectivity-at-age for fleet $`f`$,
- $`\text{Sel}_{p,r,y = 1,\tau,a,s,f}^{\text{Ret}}`$ is the retention
  selectivity-at-age for fleet $`f`$,
- $`\delta_{r,\tau,f}`$ is the discard mortality rate for fleet $`f`$ in
  region $`r`$ and season $`\tau`$,
- $`n_f`$ is the number of fishing fleets,
- $`\psi_{p,r,y,s}`$ describes the recruitment sex-ratio,
- $`\zeta_{p,r}`$ apportions the global recruitment parameter across
  regions (estimated using a multinomial logit transform to ensure
  proportions sum to one).

Because the equilibrium calculation above is a purely deterministic
(non-stochastic) projection, $`\text{rinit}_p`$ is treated as the median
of the assumed lognormal recruitment process (consistent with how
$`\mu_p^{\text{Rec}}`$ is interpreted elsewhere; see Recruitment
Processes), and the same lognormal bias-correction term used for
recruitment deviations is applied here as a static offset rather than
about an estimated deviation. This keeps the equilibrium age structure
on a scale consistent with the rest of the recruitment process even
though no annual deviation is estimated at initialization. The same
correction is applied when the operating model constructs an equivalent
equilibrium seed during closed-loop simulation
([`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)’s
`rinit_input` pathway), so fitted and simulated equilibria remain on a
consistent scale.

The plus group ($`a_{+}`$) of the initial population is then computed
as:

``` math
N_{p,r,a_{+},s}^{'} = N_{p,r,a_{+} - 1,s}^{'}\dfrac{\exp\left( - Z_{p,r,a = a_{+} - 1,s}^{'} \right)}{1 - \exp\left( - Z_{p,r,a = a_{+},s}^{'} \right)}
```

However, this scalar geometric series solution assumes that the plus
group accumulates in a closed system. Therefore, when movement dynamics
are present, this solution does not correctly accumulate individuals
into the plus group.

To address this, additional methods are provided to explicitly
incorporate movement dynamics into the plus group calculation. In
particular, the initial population can be derived by iterating the
population to equilibrium. An exponential decay model is used to first
initialize the age structure at the first iteration:

``` math
\begin{matrix}
N_{p,r,a,s}^{'} = \left\{ \begin{matrix}
\mu_p^{\text{RecInit}}\,\psi_{p,r,y = 1,s}\,\zeta_{p,r}, & \text{if }a = 1 \\
\mu_p^{\text{RecInit}}\,\psi_{p,r,y = 1,s}\,\zeta_{p,r}\,\exp\left( - \sum_{a = 1}^{n_{a}}Z_{p,r,a,s}^{'} \right), & \text{if }a > 1 \\
\end{matrix} \right.\  \\
\end{matrix}
```

The initialized age structure is then iterated forward to equilibrium by
applying recruitment, followed by the seasonal transition (movement and
mortality) and ageing. That transition is the same operator
$`\mathbf{\Phi}`$ used by the main projection, so equilibrium
initialization inherits whatever `move_timing` sequencing the model is
configured with (see Population Projection section for equations).

While the iterative method correctly accumulates the plus group when
movement is present, it can be computationally inefficient. Therefore,
`SPoRC` enables users to compute the plus group using the matrix
formulation of the geometric series, which correctly accounts for
movement processes. Let $`\mathbf{\Phi}_{p,\tau,a,s}^{'}`$ denote the
seasonal transition operator of the Population Projection section
evaluated at the equilibrium mortality $`\mathbf{Z}_{p,a,s}^{'}`$, and
define the year-long transition for age $`a`$ by composing seasons in
order (season 1 applied first, hence rightmost):

``` math
\mathbf{G}_{p,s}^{a} = \left( \mathbf{\Phi}_{p,n_\tau,a,s}^{'} \right)^{T}\left( \mathbf{\Phi}_{p,n_\tau - 1,a,s}^{'} \right)^{T}\cdots\left( \mathbf{\Phi}_{p,1,a,s}^{'} \right)^{T}
```

The transposes appear because $`\mathbf{\Phi}`$ is stored in the
row-vector convention $`\left(\mathbf{N}\right)^{T}\mathbf{\Phi}`$ used
above, whereas the recursion below is written as a column operator
acting on $`\mathbf{N}`$. The population is projected forward to the
penultimate age ($`\mathbf{N}_{p,a_{+} - 1,s}^{'}`$), and the
penultimate age is then projected forward through one more year:

``` math
\mathbf{X}_{p,s} = \mathbf{G}_{p,s}^{a_{+} - 1}\,\mathbf{N}_{p,a_{+} - 1,s}^{'}
```

where $`\mathbf{X}_{p,s}`$ represents the culmination of processes
applied to the penultimate age. $`\mathbf{G}_{p,s}^{a_{+}}`$ plays the
role of the survival ratio in the scalar series, representing the
combined effects of survival and movement on the plus group over a full
year, so the plus group solution incorporating movement is:

``` math
\mathbf{N}_{p,a_{+},s}^{'} = \left( \mathbf{I -}\mathbf{G}_{p,s}^{a_{+}} \right)^{- 1}\mathbf{X}_{p,s}
```

Under $`\text{move\_timing} = 0`$ and a single season,
$`\mathbf{G}_{p,s}^{a}`$ is
$`\text{diag}\left( \exp\left( - \mathbf{Z}_{p,a,s}^{'} \right) \right)\left( \mathbf{M}_{p,y = 1,a,s} \right)^{T}`$.
When only a single region is modeled or no movement occurs (i.e., an
identity matrix), the matrix formulation simplifies to the standard
scalar geometric series solution under every `move_timing`.

Following the definition of equilibrium age structure, initial age
deviations can be applied:

``` math
\begin{matrix}
N_{p,r,y = 1,a \neq 1,s} = N'_{p,r,a \neq 1,s}\text{exp}\left( \epsilon_{p,r,i}^{\text{Init}} \right) \\
\end{matrix}
```

where $`N_{p,r,y = 1,a \neq 1,s}`$ represents the numbers-at-age in the
first model year and season except for recruits ($`a \neq 1`$). These
values can be treated as a stochastic process by applying multiplicative
lognormal deviations $`\epsilon_{p,r,i}^{\text{Init}}`$ to the initial
equilibrium age structure. Note that the index $`i`$ is introduced
because users can determine whether initial age deviations are estimated
up to the penultimate age class, or across all classes including the
plus group.

Finally, the initial age structure can be specified as fully free
(`init_age_strc = "free"`), in which case no equilibrium is projected at
all. The deviations are then the initial numbers themselves rather than
multipliers on an equilibrium:

``` math
N_{p,r,y = 1,a \neq 1,s} = \psi_{p,r,y=1,s}\exp\left( \epsilon_{p,r,i}^{\text{Init}} \right)
```

with the first age still taken from the recruitment process. Because the
deviations are on the scale of log-numbers rather than log-ratios about
an equilibrium, the initial condition carries no information about
$`\mu_p^{\text{RecInit}}`$ or the assumed initial fishing mortality; any
penalty applied to $`\epsilon^{\text{Init}}`$ under this option acts as
a prior on log initial abundance directly. This matches the convention
of assessments in which initial numbers-at-age are freely estimated
parameters rather than deviations from an equilibrium.

The initial age deviations carry a sex dimension governed by
`InitDevs_sex_spec`. Under `"est_shared_s"` (the default) one age curve
$`\epsilon_{p,r,i}^{\text{Init}}`$ is read by every sex. Under
`"est_all"` each sex carries its own curve
$`\epsilon_{p,r,i,s}^{\text{Init}}`$, which under the free option gives

``` math
N_{p,r,y = 1,a \neq 1,s} = \psi_{p,r,y=1,s}\exp\left( \epsilon_{p,r,i,s}^{\text{Init}} \right)
```

#### Recruitment Processes

`SPoRC` holds recruitment either about a mean parameter or on a
stock-recruit curve. Recruitment can be specified to arise about a mean
parameter ($`\mu_p^{\text{Rec}}`$):

``` math
\begin{matrix}
N_{p,r,y,\tau = 1,a = 1,s} = \mu_p^{\text{Rec}}\exp\left( \epsilon_{p,r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)\chi_{p,\tau = 1}\psi_{p,r,y,s}\zeta_{p,r} \\
\end{matrix}
```

where $`\epsilon_{p,r,y}`$ are annual, lognormally distributed
recruitment deviations with a lognormal bias correction term
($`\frac{\sigma_{\text{Rec}}^{2}}{2}b_{y}`$), with $`b_{y}`$
representing the bias correction ramp from Methot and Taylor (2011), and
$`\chi_{p,\tau}`$ is the proportion of annual recruitment assigned to
season $`\tau`$ for population $`p`$ (with
$`\sum_\tau \chi_{p,\tau} = 1`$). For seasons $`\tau > 1`$, recruits are
added to the existing numbers at age 1:

``` math
N_{p,r,y,\tau > 1,a = 1,s} = N_{p,r,y,\tau > 1,a = 1,s} + \text{TotalRec}_{p,r,y} \cdot \chi_{p,\tau} \cdot \psi_{p,r,y,s}
```
where $`\text{TotalRec}_{p,r,y}`$ is the total annual recruitment
(before seasonal apportionment) for population $`p`$ in region $`r`$ and
year $`y`$.

Recruitment can also be specified to arise from a Beverton-Holt stock
recruitment function to invoke density-dependent population dynamics,
following the steepness parameterization (Mace and Doonan, 1988).
Localized density-dependent recruitment is defined as:

``` math
\begin{matrix}
N_{p,r,y,\tau = 1,a = 1,s} = \dfrac{4\mu_p^{\text{Rec}}\zeta_{p,r}h_{p,r}{\text{effSSB}}_{p,y - RecLag}}{\left( 1 - h_{p,r} \right)\text{SSB0}_{p,r} + \left( 5h_{p,r} - 1 \right){\text{effSSB}}_{p,y - RecLag}}\exp\left( \epsilon_{p,r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)\chi_{p,\tau=1}\psi_{p,r,y,s} \\
\end{matrix}
```

while global density-dependent recruitment can be defined as:

``` math
\begin{matrix}
N_{p,r,y,\tau = 1,a = 1,s} = \dfrac{4\mu_p^{\text{Rec}}\zeta_{p,r}h_p\sum_{r}^{}{SSB}_{p,r,y - RecLag}}{(1 - h_p)\sum_{r}^{}{SSB0_{p,r}} + (5h_p - 1)\sum_{r}^{}{SSB}_{p,r,y - RecLag}}\exp\left( \epsilon_{p,r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)\chi_{p,\tau=1}\psi_{p,r,y,s} \\
\end{matrix}
```

where $`\mu_p^{\text{Rec}}`$ under this parameterization is the virgin
unfished recruitment for population $`p`$, $`h_{p,r}`$ (or $`h_p`$) is
the steepness parameter representing the fraction of
$`\mu_p^{\text{Rec}}\zeta_{p,r}`$ that would be produced when at 20% of
$`\text{SSB0}_{p,r}`$ (or $`\sum_{r}^{}{SSB0_{p,r}}`$). The steepness
parameter is constrained to be between values of 0.2 and 1 and are
estimated in bounded logit space. $`\text{SSB0}_{p,r}`$ is a derived
variable that represents the unfished spawning stock biomass.
$`{SSB}_{p,r,y - RecLag}`$ is the spawning stock biomass for population
$`p`$ in region $`r`$, and $`\text{effSSB}_{p,y}`$ is the effective
spawning stock biomass (see Spawning Biomass section below).

A Ricker stock-recruit relationship can alternatively be specified
(`rec_model = "ricker_rec"`), written in a depletion form. For localized
density dependence:

``` math
N_{p,r,y,\tau = 1,a = 1,s} = \mu_p^{\text{Rec}}\zeta_{p,r}\frac{\text{effSSB}_{p,y - RecLag}}{\text{SSB0}_{p,r}}\exp\left( \alpha_{p,r}\left\lbrack 1 - \frac{\text{effSSB}_{p,y - RecLag}}{\text{SSB0}_{p,r}} \right\rbrack \right)\exp\left( \epsilon_{p,r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)\chi_{p,\tau=1}\psi_{p,r,y,s}
```

with the global form defined analogously by replacing the regional
spawning biomass and $`\text{SSB0}`$ with their sums across regions. The
log-slope is derived from steepness as:

``` math
\alpha_{p,r} = \log\left( \frac{4h_{p,r}}{1 - h_{p,r}} \right)
```

so the curve passes through $`(\text{SSB0}, R_0)`$ by construction and
carries the same compensation ratio as a Beverton-Holt at the same
$`h`$. This is not the textbook steepness definition
$`R(0.2\,\text{SSB0}) = h R_0`$: the Ricker here yields
$`R(0.2\,\text{SSB0})/R_0 = 0.2\left(4h/(1-h)\right)^{0.8}`$, which
exceeds $`h`$ and is not bounded above by 1. Steepness values are
therefore not interchangeable between the Beverton-Holt and Ricker
forms, and a steepness prior calibrated for one should not be reused for
the other without translation.

For both stock-recruit forms, the biological inputs to unfished spawning
biomass per recruit (weight-at-age, maturity, natural mortality, and
movement), and hence to $`\text{SSB0}`$ and the scale of the curve, are
taken from a single reference year set by `SR_ref_yr`. The default is
the first model year, which is the model’s long-standing behavior;
setting `SR_ref_yr` to the terminal year conditions the curve on
terminal-year biologicals instead, a convention several assessments use.
With time-varying weight-at-age the two choices imply different
$`\text{SSB0}`$ and therefore different depletion scales.

The spawning stock biomass is the product of numbers-at-age, spawning
weight-at-age, and maturity-at-age for females in the spawning season
$`\tau^{spawn}`$:

``` math
\begin{matrix}
SSB_{p,r,y} = \sum_{a = 1}^{a_{+}}{N_{p,r,y,\tau^{spawn},a,s = 1}^{spawn}W_{p,r,y,\tau^{spawn},a,s = 1}^{spawn}\text{Mat}_{p,r,y,\tau^{spawn},a,s = 1}} \\
\end{matrix}
```

where $`N^{spawn}`$ is the numbers-at-age propagated from the start of
the spawning season to the spawning point $`t^{spawn}`$ within it. That
propagation depends on how movement and mortality are sequenced and is
defined in the Spawning Biomass Timing section below. Under the default
$`\text{move\_timing} = 0`$ it is simply
$`N_{p,r,y,\tau^{spawn},a,s = 1}\exp\left(-Z_{p,r,y,\tau^{spawn},a,s=1} \cdot t^{spawn}\right)`$
evaluated after that season’s movement has been applied, which is the
familiar form.

For single-sex models, SSB is multiplied by 0.5 to obtain female-only
spawning biomass.

Note that $`RecLag`$ denotes the delay (in seasons) between spawning and
when recruits enter the population, and is user-specified as
$`RecLag \geq 1`$ (the classic case above) or $`RecLag = 0`$ (age-0
recruitment, described below). For $`RecLag \geq 1`$: if
$`y \leq RecLag`$ (i.e. there is not yet enough model history to look
back $`RecLag`$ seasons), `SPoRC` utilizes the initial equilibrium
spawning biomass instead of $`{SSB}_{p,r,y - RecLag}`$ to compute
deterministic recruitment. That equilibrium is evaluated at the initial
fishing mortality, so it equals $`\text{SSB0}_{p,r}`$ whenever the model
starts unfished.

##### Stock-Recruit Curve as a Penalty

Under both density-dependent forms above the curve generates
recruitment, which makes $`\epsilon_{p,r,y}^{\text{Rec}}`$ the
stock-recruit residual by construction. A third arrangement
(`rec_model = "mean_rec"` with `sr_penalty = "bh"` or `"ricker"`) leaves
recruitment arising about the mean parameter exactly as in the first
equation of this section,

``` math
R_{p,r,y} = \mu_p^{\text{Rec}}\zeta_{p,r}\exp\left( \epsilon_{p,r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)
```

where $`R_{p,r,y} \equiv \text{TotalRec}_{p,r,y}`$ is the realized total
annual recruitment, and evaluates a stock-recruit curve alongside the
dynamics without ever advancing them. The curve supplies a prediction
$`{\widehat{R}}_{p,r,y}`$ (reported as `SR_pred`, with dimensions
population by region by year) from the Beverton-Holt or Ricker
expression given above, driven by the same
$`\text{effSSB}_{p,y - RecLag}`$ (or
$`\sum_{r}^{}{SSB}_{p,r,y - RecLag}`$ under global density dependence)
and carrying no recruitment deviation. Numbers-at-age are advanced by
$`R_{p,r,y}`$ and never by $`{\widehat{R}}_{p,r,y}`$, so the deviations
remain free about the mean rather than becoming residuals about the
curve. The curve enters the model only through the log residual

``` math
\xi_{p,r,y}^{\text{SR}} = \log R_{p,r,y} - \log{\widehat{R}}_{p,r,y}
```

which enters the objective as a penalty (see the Stock-Recruit Residual
section under Process Error Penalties). Note that
$`\xi_{p,r,y}^{\text{SR}}`$ is a derived quantity rather than an
estimated parameter, which is what separates this arrangement from the
two above: a weakly determined relationship informs the recruitment
series at a cost set by the penalty instead of dictating it, the
convention several existing assessments use. The two statements are
mutually exclusive, and a stock-recruit penalty requested alongside
`rec_model = "bh_rec"` or `"ricker_rec"` is rejected at setup, since the
curve already generates recruitment there and the residual would be
penalized twice.

The scale of the curve, and with it $`\text{SSB0}_{p,r}`$ and the
depletion at which the curve is evaluated, is set by `sr_R0_spec`:

- `"shared"` (the default) takes the scale from $`\mu_p^{\text{Rec}}`$,
  which under mean recruitment is the level of the recruitment series
  itself. A single parameter carries both the recruitment level and the
  curve, which is better posed than `"est"`, and it anchors the curve on
  mean recruitment rather than on an unfished level.
- `"est"` gives the curve its own estimated scale $`\mu_p^{\text{SR}}`$,
  identified by the curve fit alone. This reproduces templates that
  carry separate mean-recruitment and unfished-recruitment parameters.
  Nothing ties $`\mu_p^{\text{SR}}`$ to $`\mu_p^{\text{Rec}}`$, so the
  two are free to slide against one another.
- `"rinit"` takes the scale from $`\mu_p^{\text{RecInit}}`$, the initial
  equilibrium recruitment of the Population Initialization section, so
  one parameter sets both the unfished age structure and the curve. This
  requires a separately estimated initialization recruitment
  (`use_rinit = 1`) and is the usual ADMB arrangement.

##### Age-0 Recruitment ($`RecLag = 0`$)

When $`RecLag = 0`$, recruitment for year $`y`$ is driven by that same
year’s own spawning biomass rather than a prior year’s. The
Beverton-Holt and Ricker equations above still apply, but with
$`\text{effSSB}_{p,y}`$ (or $`\sum_r SSB_{p,r,y}`$ for global density
dependence) in place of $`\text{effSSB}_{p,y - RecLag}`$, and recruits
enter starting at the spawning season $`\tau^{spawn}`$ rather than
season 1, with any remaining seasonal share ($`\tau > \tau^{spawn}`$)
added to the existing numbers at age 1 exactly as in the $`\tau > 1`$
equation above.

Because $`\text{effSSB}_{p,y}`$/$`SSB_{p,r,y}`$ is not knowable until
season $`\tau^{spawn}`$ is actually reached within year $`y`$, this
timing constraint is enforced structurally rather than left to the user:

- $`\chi_{p,\tau} = 0`$ for every season before $`\tau^{spawn}`$
  (validated at setup; recruits cannot predate the spawning event that
  produced them), and
- $`\text{Mat}_{p,r,y,\tau,a = 1,s} = 0`$ for all $`\tau`$ (validated at
  setup; age-0 fish cannot be mature), which guarantees the $`a = 1`$
  term in the $`SSB_{p,r,y}`$ sum is always zero regardless of whether
  this year’s recruits have been added to
  $`N_{p,r,y,\tau^{spawn},a=1,s}`$ yet at the point $`SSB_{p,r,y}`$ is
  evaluated.

$`\text{SSB0}_{p,r}`$ itself is a pure per-recruit, equilibrium quantity
and does not depend on $`RecLag`$ which is the same value is used
whether recruitment is lagged or age-0. Unlike the $`RecLag \geq 1`$
case, there is no burn-in substitution of $`\text{SSB0}_{p,r}`$ for
early years: since $`RecLag = 0`$, $`SSB_{p,r,y}`$ (this year’s own
survivor biomass) is always available by the time it is needed,
including in year 1.

Because age-0 recruits are inserted in the middle of the spawning season
rather than at its start, `move_timing` also determines whether they
need a movement step of their own. Under $`\text{move\_timing} = 0`$ the
season’s movement has already been applied by the time $`SSB_{p,r,y}`$
is evaluated, so recruits inserted immediately afterwards would
otherwise miss it entirely; when $`\text{do\_recruits\_move} = 1`$ they
are therefore given that season’s movement step explicitly, and they
then experience the rest of the season’s mortality like any other
seasonal recruit pulse. Under $`\text{move\_timing}`$ 1 and 2 movement
has not yet occurred at that point in the season, so the newly inserted
recruits are carried by the end-of-season transition $`\mathbf{\Phi}`$
along with every other age and no separate catch-up step is applied.

##### Effective Spawning Biomass and Multi-Population Dynamics

When multiple populations are modeled ($`n_p > 1`$), effective spawning
biomass at each population’s natal region accounts for stray
contributions from other populations:

``` math
\begin{matrix}
\text{effSSB}_{p,y} = SSB_{p, r^{\text{natal}}_p, y} + \sum_{p' \neq p} \frac{\phi_{p',y}}{npop_r} \cdot SSB_{p', r^{\text{natal}}_p, y}
\end{matrix}
```

where $`r^{\text{natal}}_p`$ is the natal region of population $`p`$,
$`\phi_{p',y}`$ is the stray rate of population $`p'`$ (the fraction of
its spawning biomass contributing to non-natal regions), and the sum is
taken over all other populations $`p' \neq p`$. For a single population,
$`\text{effSSB}_{1,y} = \sum_r SSB_{1,r,y}`$. Note that $`npop_r`$ is
the number of populations in a given region, where the contribution of
$`\phi_{p',y}`$ is split evenly among populations.

##### Single-Season Spawning Movement

When $`n_\tau = 1`$ and $`n_p > 1`$, a separate spawning movement matrix
$`\mathbf{M}^{spawn}_{p,y,a,s}`$ is applied to both fished and unfished
numbers-at-age prior to computing spawning biomass quantities,
representing natal homing of individuals to their spawning grounds. It
composes with, rather than replaces, whatever within-season movement
`move_timing` implies: it acts on the numbers already propagated to the
spawning point,

``` math
\left(\mathbf{N}^{homed}_{p,y,a,s}\right)^T = \left(\mathbf{N}^{spawn}_{p,y,\tau^{spawn},a,s}\right)^T \mathbf{M}^{spawn}_{p,y,a,s}
```

with $`\mathbf{N}^{spawn}`$ as defined in the Spawning Biomass Timing
section below. The ordering relative to the $`t^{spawn}`$ mortality
discount differs by timing, and the two do not commute once movement
redistributes fish: under $`\text{move\_timing} = 0`$ natal homing is
applied to the post-movement numbers and the $`t^{spawn}`$ discount is
taken afterwards, at the destination region, whereas under
$`\text{move\_timing}`$ 1 and 2 the propagation to the spawning point
already carries that discount and natal homing is applied to the result.

This additional movement is applied only for spawning biomass
calculations and does not alter the numbers-at-age array used for
subsequent mortality and movement processes.

#### Population Projection

Following recruitment processes, the population is projected forward. In
the spatial model, each season advances the population under both
Markovian movement and total mortality. Movement is described by a
first-order Markov matrix $`\mathbf{M}_{p,y,\tau,a,s}`$ acting on the
numbers-at-region vector,
$`\left( \mathbf{N}_{p,y,\tau,a,s} \right)^{T}\mathbf{M}_{p,y,\tau,a,s}`$.
In a single-region case, no movement is applied (i.e.,
$`\mathbf{M}_{p,y,\tau,a,s}`$ is an implied identity matrix). For each
population, year, season, age, and sex combination, the movement matrix
specifies bulk-transfer coefficients. Where within the season that
transfer happens, before mortality, after it, or continuously alongside
it, is set by `move_timing`, and is what the rest of this section makes
explicit.

Movement and mortality are therefore combined into a single seasonal
transition operator $`\mathbf{\Phi}_{p,y,\tau,a,s}`$. Writing
$`\mathbf{Z}_{p,y,\tau,a,s}`$ for the vector of seasonal total mortality
across regions,
$`\mathbf{s}_{p,y,\tau,a,s} = \exp\left( -\mathbf{Z}_{p,y,\tau,a,s} \right)`$
for seasonal survival, and $`\dot{\mathbf{Q}}_{p,y,\tau,a,s}`$ for the
CTMC generator:

``` math
\mathbf{\Phi}_{p,y,\tau,a,s} = \begin{cases}
\mathbf{M}_{p,y,\tau,a,s}\,\text{diag}\left( \mathbf{s}_{p,y,\tau,a,s} \right) & \text{move\_timing} = 0 \\[4pt]
\text{diag}\left( \mathbf{s}_{p,y,\tau,a,s} \right)\mathbf{M}_{p,y,\tau,a,s} & \text{move\_timing} = 1 \\[4pt]
\left[ \exp\left( \mathbf{\Lambda}_{p,y,\tau,a,s} \right) \right]^{T} & \text{move\_timing} = 2 \\
\end{cases}
```

where

``` math
\mathbf{\Lambda}_{p,y,\tau,a,s} = \dot{\mathbf{Q}}_{p,y,\tau,a,s}^{T}\,\Delta\tau - \text{diag}\left( \mathbf{Z}_{p,y,\tau,a,s} \right)
```

is the combined movement-mortality generator for the season, expressed
in the column convention ($`\dot{\mathbf{Q}}`$ is stored row-wise, hence
the transpose). Note that $`\mathbf{\Lambda}`$ already carries the
season duration in both of its terms ($`\Delta\tau`$ explicitly on the
generator, and implicitly in $`\mathbf{Z}`$, which is itself a seasonal
rate), so the season is parameterized on the unit interval throughout:
propagating a fraction $`t`$ of the way through the season means
$`\exp\left( \mathbf{\Lambda}\,t \right)`$ with $`t \in [0,1]`$, and a
full season is $`t = 1`$.

Under $`\text{move\_timing} = 0`$ (sequential; movement then mortality)
individuals move at the start of the season and then experience
mortality in the destination region. Under $`\text{move\_timing} = 1`$
(sequential; mortality then movement) they experience mortality in the
origin region and move at the end of the season. Under
$`\text{move\_timing} = 2`$ (continuous) movement and mortality act
simultaneously, with the generator and mortality rate combined inside a
single matrix exponential; $`\mathbf{M}`$ does not appear at all, and
this option requires an estimated CTMC generator
($`\text{move\_type} = 1`$). The three coincide when $`\mathbf{Z}`$ is
constant across regions (a scalar multiple of the identity commutes with
the generator), when $`\mathbf{Z} = \mathbf{0}`$, and when
$`\mathbf{M} = \mathbf{I}`$ (equivalently
$`\dot{\mathbf{Q}} = \mathbf{0}`$).

For seasons within a year ($`\tau < n_\tau`$), individuals advance to
the next season at the same age:

``` math
\left( \mathbf{N}_{p,y,\tau + 1,a,s} \right)^{T} = \left( \mathbf{N}_{p,y,\tau,a,s} \right)^{T}\mathbf{\Phi}_{p,y,\tau,a,s}
```

At the end of the final season ($`\tau = n_\tau`$), individuals advance
in age:

``` math
\left( \mathbf{N}_{p,y + 1,1,a + 1,s} \right)^{T} = \left( \mathbf{N}_{p,y,n_\tau,a,s} \right)^{T}\mathbf{\Phi}_{p,y,n_\tau,a,s},\quad\text{for }1 \leq a < a_{+}
```

``` math
\left( \mathbf{N}_{p,y + 1,1,a_{+},s} \right)^{T} = \left( \mathbf{N}_{p,y + 1,1,a_{+},s} \right)^{T} + \left( \mathbf{N}_{p,y,n_\tau,a_{+},s} \right)^{T}\mathbf{\Phi}_{p,y,n_\tau,a_{+},s}
```

When $`\text{move\_timing} = 0`$ these reduce to the elementwise
sequential form, with movement applied first as
$`\mathbf{N}_{p,y,\tau,a,s} = \left( \mathbf{N}_{p,y,\tau,a,s} \right)^{T}\mathbf{M}_{p,y,\tau,a,s}`$
and mortality applied afterwards as:

``` math
N_{p,r,y,\tau + 1,a,s} = N_{p,r,y,\tau,a,s}\exp\left( - Z_{p,r,y,\tau,a,s} \right)
```

``` math
N_{p,r,y + 1,1,a + 1,s} = N_{p,r,y,n_\tau,a,s}\exp\left( - Z_{p,r,y,n_\tau,a,s} \right),\quad\text{for }1 \leq a < a_{+}
```

``` math
N_{p,r,y + 1,1,a_{+},s} = N_{p,r,y+1,1,a_{+},s} + N_{p,r,y,n_\tau,a_{+},s}\exp\left( - Z_{p,r,y,n_\tau,a_{+},s} \right)
```

If recruits do not move ($`\text{do\_recruits\_move} = 0`$),
$`\mathbf{M}_{p,y,\tau,a=1,s} = \mathbf{I}`$ and
$`\dot{\mathbf{Q}}_{p,y,\tau,a=1,s} = \mathbf{0}`$, so
$`\mathbf{\Phi}_{p,y,\tau,a=1,s}`$ reduces to
$`\text{diag}\left( \mathbf{s}_{p,y,\tau,a=1,s} \right)`$ under all
three cases.

##### Where in the Season Each Quantity Is Evaluated

Because $`\mathbf{\Phi}`$ resolves movement at a different point in the
season under each option, every quantity that is observed partway
through a season has to be evaluated consistently. The table below
summarizes where each one is taken; the equations follow in the sections
indicated.

| Quantity | $`\text{move\_timing} = 0`$ | $`\text{move\_timing} = 1`$ | $`\text{move\_timing} = 2`$ |
|----|----|----|----|
| Spawners, $`\mathbf{N}^{spawn}`$ | Post-movement location, discounted by $`t^{spawn}`$ mortality | Pre-movement (origin) location, discounted by $`t^{spawn}`$ mortality | Partially redistributed: $`\exp(\mathbf{\Lambda}\,t^{spawn})`$ |
| Catch and discards | Baranov at post-movement location | Baranov at pre-movement location | Season-integrated abundance $`\mathbf{N}^{\int}`$ (spatial Baranov) |
| Fishery index | Post-movement $`\mathbf{N}`$, no discount | Pre-movement $`\mathbf{N}`$, no discount | $`\mathbf{N}^{\int}`$ |
| Survey index | Post-movement $`\mathbf{N}`$, discounted by $`t^{srv}`$ mortality | Pre-movement $`\mathbf{N}`$, discounted by $`t^{srv}`$ mortality | Snapshot $`\exp(\mathbf{\Lambda}\,t^{srv})\mathbf{N}`$ |
| Tag recaptures | Baranov on the tag cohort | Baranov on the tag cohort | Season-integrated tag abundance |

Under continuous movement, note that the survey index is treated
differently from catch and the fishery index. A survey is a snapshot at
an instant within the season, so it uses the partial propagation
$`\exp(\mathbf{\Lambda}\,t^{srv})`$, whereas catch and the fishery index
accumulate over the whole season, so they use the integral of that
propagation. Under the sequential timings the distinction does not
arise, because movement has already been resolved and the fish are
stationary for the rest of the season.

##### Spawning Biomass Timing

Spawning biomass is computed from the population propagated a fraction
$`t^{spawn}`$ into the spawning season, consistently with the
sequencing:

``` math
\left( \mathbf{N}_{p,y,\tau^{spawn},a,s}^{spawn} \right)^{T} = \begin{cases}
\left( \left( \mathbf{N}_{p,y,\tau^{spawn},a,s} \right)^{T}\mathbf{M}_{p,y,\tau^{spawn},a,s} \right)\text{diag}\left( \exp\left( - t^{spawn}\mathbf{Z}_{p,y,\tau^{spawn},a,s} \right) \right) & \text{move\_timing} = 0 \\[4pt]
\left( \mathbf{N}_{p,y,\tau^{spawn},a,s} \right)^{T}\text{diag}\left( \exp\left( - t^{spawn}\mathbf{Z}_{p,y,\tau^{spawn},a,s} \right) \right) & \text{move\_timing} = 1 \\[4pt]
\left( \mathbf{N}_{p,y,\tau^{spawn},a,s} \right)^{T}\left[ \exp\left( \mathbf{\Lambda}_{p,y,\tau^{spawn},a,s}\,t^{spawn} \right) \right]^{T} & \text{move\_timing} = 2 \\
\end{cases}
```

so spawners are at their post-movement locations under
$`\text{move\_timing} = 0`$, at their pre-movement locations under
$`\text{move\_timing} = 1`$, and partially redistributed under
$`\text{move\_timing} = 2`$. No extra convention is imposed to make the
three agree: each timing has its own natural spawning state, and the
$`\text{move\_timing} = 0`$ case reproduces the historical `SPoRC`
calculation exactly. The same $`\mathbf{N}^{spawn}`$ is used for
$`SSB_{p,r,y}`$, total biomass, and the unfished
$`\text{Dynamic\_SSB0}_{p,r,y}`$ (the last with $`\mathbf{Z}`$ replaced
by natural mortality alone).

##### Catch Under Continuous Movement

The Baranov catch equation assumes individuals remain in one region for
the whole season, which does not hold under $`\text{move\_timing} = 2`$.
Catch is instead taken from the season-integrated abundance:

``` math
\mathbf{N}_{p,y,\tau,a,s}^{\int} = \left[ \int_{0}^{1}\exp\left( \mathbf{\Lambda}_{p,y,\tau,a,s}\,\upsilon \right)d\upsilon \right]\mathbf{N}_{p,y,\tau,a,s}
```

``` math
C_{p,r,y,\tau,a,s,f} = \text{retFmort}_{p,r,y,\tau,a,s,f} \cdot N_{p,r,y,\tau,a,s}^{\int}
```

``` math
D_{p,r,y,\tau,a,s,f} = \text{discFmort}_{p,r,y,\tau,a,s,f} \cdot N_{p,r,y,\tau,a,s}^{\int}
```

The integral runs over the unit interval rather than over
$`[0, \Delta\tau]`$ because, as noted above, $`\mathbf{\Lambda}`$
already carries the season duration; $`\upsilon`$ is elapsed fraction of
the season, not elapsed time.

When $`\dot{\mathbf{Q}} = \mathbf{0}`$ the regions decouple and this
reduces to
$`N_{p,r,y,\tau,a,s}\left( 1 - \exp\left( - Z_{p,r,y,\tau,a,s} \right) \right)/Z_{p,r,y,\tau,a,s}`$,
recovering the standard Baranov form used under
$`\text{move\_timing} \in \{0,1\}`$, which makes explicit that
$`\mathbf{N}^{\int}`$ is simply what the familiar $`(1 - e^{-Z})/Z`$
factor was always computing, namely the abundance accumulated over the
season. The same season-integrated abundance is used for the fishery
index and for predicted tag recaptures when movement is specified to be
continuous.

$`Z_{p,r,y,\tau,a,s}`$ denotes the seasonal total instantaneous
mortality rate and is defined as the combination of natural mortality
($`\text{Natmort}_{p,r,y,a,s}`$) scaled by seasonal duration
$`\Delta\tau`$, retained fishing mortality
($`\text{retFmort}_{p,r,y,\tau,a,s,f}`$), and dead discard fishing
mortality ($`\text{discFmort}_{p,r,y,\tau,a,s,f}`$):

``` math
\begin{matrix}
Z_{p,r,y,\tau,a,s} = \text{Natmort}_{p,r,y,a,s} \cdot \Delta\tau + \sum_{f}^{}\left[\text{retFmort}_{p,r,y,\tau,a,s,f} + \text{discFmort}_{p,r,y,\tau,a,s,f}\right] \\
\end{matrix}
```

where the retained and dead discard fishing mortality rates at age are:

``` math
\text{retFmort}_{p,r,y,\tau,a,s,f} = \text{Fmort}_{r,y,\tau,f} \cdot \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}} \cdot \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}}
```

``` math
\text{discFmort}_{p,r,y,\tau,a,s,f} = \text{Fmort}_{r,y,\tau,f} \cdot \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}} \cdot \left(1 - \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}}\right) \cdot \delta_{r,y,\tau,f}
```

Here, $`\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}}`$ is the total
fishery selectivity (governing encounter probability),
$`\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}}`$ is the retention
selectivity (governing the probability of retention given encounter),
and $`\delta_{r,y,\tau,f}`$ is the discard mortality rate for fleet
$`f`$. Only the dead fraction of discards contributes to total
mortality. The seasonal instantaneous fishing mortality rate is defined
as:

``` math
\begin{matrix}
\text{Fmort}_{r,y,\tau,f} = \mu_{r,\tau,f}^{\text{Fsh}}\text{exp}\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} \right) \\
\end{matrix}
```

where $`\text{Fmort}_{r,y,\tau,f}`$ is parameterized based on lognormal
deviations ($`\epsilon_{r,y,\tau,f}^{\text{Fsh}})`$ about a mean fishing
mortality parameter for a given region, season, and fishery fleet
($`\mu_{r,\tau,f}^{\text{Fsh}}`$). When no catch data are available for
a given region, season, and fleet, fishing mortality is set to zero.

Under the default `ln_F_mean_spec = "est"` the mean is estimated and the
deviations are departures from it. `ln_F_mean_spec = "fix"` instead
fixes the mean at its starting value (zero on the log scale unless
supplied), so that

``` math
\text{Fmort}_{r,y,\tau,f} = \exp\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} \right)
```

and the deviations are annual log fishing mortality outright, estimated
as free parameters. Because the deviation penalty is then the only
statement about the level of $`F`$, this parameterization must be paired
with a penalty that leaves the level free:
`Fdev_pen_center = "own_mean"` (deviations penalized about their own
mean), `Fdev_model = "rw"` (only increments penalized), or
`Use_F_pen = 0`. An `"iid"` or `"ar1"` penalty centered on the fixed
zero mean would shrink the deviations toward $`F = 1`$, so that
combination is rejected at setup.

The discard mortality rate is parameterized analogously via logistic
deviations about a mean logit-scale discard mortality rate:

``` math
\delta_{r,y,\tau,f} = \text{logistic}\left(\mu_{r,\tau,f}^{\delta} + \epsilon_{r,y,\tau,f}^{\delta}\right)
```

where $`\mu_{r,\tau,f}^{\delta}`$ is the logit-scale mean discard
mortality rate and $`\epsilon_{r,y,\tau,f}^{\delta}`$ are annual
deviations. The discard mortality rate is bounded between 0 and 1.

#### Movement Processes

Movement processes can be parameterized as either an unstructured Markov
process (discrete-time Markov) or as a Continuous-time Markov chain
(CTMC) process. Movement parameterized as an unstructured Markov process
is estimated using a multinomial logit link function, with
$`n_r \times (n_r-1)`$ free parameters per stratum:

``` math
M_{p,r,k,y,\tau,a,s} = \frac{\exp(\omega_{p,r,k,y,\tau,a,s})}{\sum_k\exp(\omega_{p,r,k,y,\tau,a,s})}
```

where the reference region $`k = 1`$ is set as
$`\omega_{p,r,k=1,y,\tau,a,s} = 0`$. Under this parameterization,
movement fractions can be estimated independently for each stratum
$`(p,y,\tau,a,s)`$, or grouped into blocks to reduce the number of
parameters.

Alternatively, movement can be specified as a continuous-time Markov
chain (CTMC) process, which decomposes into diffusive and taxis
components. Diffusive processes represent undirected movement of
individuals, while taxis processes represent directed movement toward
more preferred habitat. This CTMC movement parameterization is governed
by an adjacency matrix ($`A_{r,k}`$), which defines neighboring regions
that can receive individuals within a given time step. Diffusive
processes are given by:

``` math
  \dot{D}_{p,r,k,y,\tau,a,s} =
  \begin{cases} 
  \dfrac{e^{2\theta}}{V_r}, & \text{if } A_{r,k} = 1 \text{ and } r \neq k, \\[0.5em]
  - \sum_{j \neq r} \dot{D}_{p,r,j,y,\tau,a,s}, & \text{if } r = k, \\[0.5em]
  0, & \text{otherwise.}
  \end{cases}
```

where $`\dot{D}_{p,r,k,y,\tau,a,s}`$ represents diffusion, $`\theta`$ is
the log diffusion rate, $`V_r`$ scales the log diffusion rate, such that
smaller regions have higher diffusion rates, and $`j`$ in the second
equation indexes all destinations except the source to ensure that the
rows of the matrix sum to 0, thereby conserving abundance. Taxis
processes (preference) can then be written as:

``` math
\dot{P}_{r,k,y,a,s} =
\begin{cases} 
h_{k,y,a,s} - h_{r,y,a,s}, & \text{if } A_{r,k} = 1 \text{ and } r \neq k, \\[0.5em]
- \sum_{j \neq r} \dot{P}_{r,j,y,a,s}, & \text{if } r = k, \\[0.5em]
0, & \text{otherwise.}
\end{cases}
```

Here, $`\dot{P}_{p,r,k,y,\tau,a,s}`$ represents the taxis (preference)
component of movement. Equation 1 determines local differences in the
habitat preference function, $`\mathbf{h}_{y,a,s}`$, while equation 2
ensures that the rows of the matrix sum to 0, conserving abundance.

Habitat preference can be defined flexibly as a combination of linear
effects and basis splines:

``` math
h_{p,r,y,\tau,a,s} = \sum_{k=1}^{n_k} \beta_{r,k} W_{p,r,k,y,\tau,a,s}
```

where $`\beta_{r,k}`$ are the estimated effects (incorporating linear or
spline effects), and $`W_{p,r,k,y,\tau,a,s}`$ is the design matrix.

Diffusive and taxis processes can then be combined to construct a
generator matrix:

``` math
\dot{Q}_{p,r,k,y,\tau,a,s} = \dot{D}_{p,r,k,y,\tau,a,s} + \dot{P}_{p,r,k,y,\tau,a,s}
```

Here, $`\dot{Q}_{p,r,k,y,\tau,a,s}`$ represents instantaneous movement
rates. The generator matrix $`\dot{\mathbf{Q}}_{p,y,\tau,a,s}`$ is
Metzler, such that all off-diagonal elements satisfy:

``` math
\dot{Q}_{p,y,\tau,a,s} \ge 0 \quad \text{for } r \neq k.
```

Under the sequential timings ($`\text{move\_timing}`$ 0 and 1), the
instantaneous movement matrix is converted to movement fractions using
the matrix exponential:

``` math
\mathbf{M}_{p,y,\tau,a,s} = \exp\Big( \dot{\mathbf{Q}}_{p,y,\tau,a,s} \, \Delta t \Big)
```

where $`\Delta t`$ is the duration of the movement interval. When
`ctmc_scale_by_seasdur = 1` (default), $`\dot{\mathbf{Q}}`$ is an annual
rate and $`\Delta t = \Delta\tau`$, so the generator and $`Z`$ share
time units; when `0`, $`\Delta t = 1`$ irrespective of season duration.

Every exponential above, here and in the $`\text{move\_timing} = 2`$
operators, is evaluated according to `move_expm_nsub`. The default `0`
takes it exactly. A value a power of two $`n \ge 1`$ instead applies
$`n`$ implicit (backward Euler) steps,

``` math
\exp(\mathbf{A}) \;\approx\; \left( \mathbf{I} - \frac{\mathbf{A}}{n} \right)^{-n},
```

which is the Pade(0,1) approximant and is first order in $`1/n`$. It
exists because the reverse-mode derivative of a linear solve is another
solve, while that of a matrix exponential is much more costly, so the
gradient is several times cheaper. Because $`\mathbf{I} - \mathbf{A}/n`$
is a non-singular M-matrix its inverse is non-negative, and with
$`\mathbf{Z} = \mathbf{0}`$ the columns sum to one exactly, so movement
fractions remain a valid transition matrix at any $`n`$. The
approximation is not free, though: at $`n = 1`$ survival is
$`1/(1 + Z)`$ rather than $`e^{-Z}`$, and simulation testing puts the
resulting bias in the estimated diffusion rate at roughly 12%, biased
high.

Under $`\text{move\_timing} = 2`$ this exponential is never taken on its
own. The generator enters the seasonal operator together with mortality,
as
$`\exp\left( \dot{\mathbf{Q}}^{T}\Delta\tau - \text{diag}(\mathbf{Z}) \right)`$,
which does not factor into a movement matrix times a survival matrix
unless $`\mathbf{Z}`$ is constant across regions. $`\mathbf{M}`$ is
still computed and reported, but as a diagnostic, the movement fractions
that would apply in the absence of spatially varying mortality, rather
than as a term in the dynamics. `ctmc_scale_by_seasdur` is forced to `1`
in this case, since combining an unscaled generator with seasdur-scaled
mortality inside a single exponential is dimensionally inconsistent.

Continuous movement also requires an estimated generator: it is
available only for $`\text{move\_type} = 1`$ with
$`\text{use\_fixed\_movement} = 0`$. Unstructured multinomial-logit
movement supplies transition fractions with no guaranteed real generator
behind them (the Markov embedding problem), so `SPoRC` rejects the
combination rather than attempting a matrix logarithm.

### Observation Equations

#### Growth and the Size-Age Transition

The size-age transition $`\mathbf{A}_{p,r,y,\tau,s}^{l}`$ converts
numbers at age to numbers at length: column $`a`$ holds the probability
that a fish of age $`a`$ falls in each length bin, and sums to one. It
is either supplied as data or built inside the model from growth
parameters (`growth_model = "vb_schnute"`), in which case every fleet
gets its own, read at that fleet’s timing in the season. Four pieces go
into it: the mean length at age, the spread of length about that mean,
the treatment of the plus group, and the binning.

##### Mean length at age

Growth is von Bertalanffy in Schnute’s parameterization: instead of
$`L_{\infty}`$ and $`t_{0}`$, the curve is anchored by the mean lengths
$`L_{1}`$ and $`L_{2}`$ at two reference ages $`A_{1}`$ and $`A_{2}`$,
with $`K`$ the rate of approach to the asymptote. The parameters are
then on the scale of the data, which makes them easier to start and to
bound. Below $`A_{1}`$ the curve is not extrapolated; mean length rises
linearly from $`L_{0}`$, the lower edge of the first length bin, to
$`L_{1}`$. Writing $`x`$ for the real age (the integer age plus the
fraction of the year elapsed):

``` math
\bar{L}_{p,r,s}(x) = \begin{cases}
L_{0} + \dfrac{x}{A_{1}}\left( L_{1,p,r,s} - L_{0} \right) & x < A_{1} \\[2ex]
L_{\infty,p,r,s} + \left( L_{1,p,r,s} - L_{\infty,p,r,s} \right)\exp\left\lbrack - K_{p,r,s}\left( x - A_{1} \right) \right\rbrack & x \geq A_{1}
\end{cases}
```

where the asymptote follows from the two anchors,

``` math
L_{\infty,p,r,s} = L_{1,p,r,s} + \dfrac{L_{2,p,r,s} - L_{1,p,r,s}}{1 - \exp\left\lbrack - K_{p,r,s}\left( A_{2} - A_{1} \right) \right\rbrack}
```

Setting `growth_A2 = "Linf"` reads $`L_{2}`$ as $`L_{\infty}`$ itself.
The five parameters $`L_{1}, L_{2}, K, CV_{1}, CV_{2}`$ are estimated on
the log scale, one set per population, region and sex, or shared across
regions or sexes (`growth_spec`).

Under `growth_model = "richards"` the same curve is applied to lengths
raised to a power $`\rho`$, a sixth estimated parameter, which lets the
inflection sit anywhere rather than at the origin:

``` math
\bar{L}_{p,r,s}(x)^{\rho} = L_{\infty,p,r,s}^{\rho} + \left( L_{1,p,r,s}^{\rho} - L_{\infty,p,r,s}^{\rho} \right)\exp\left\lbrack - K_{p,r,s}\left( x - A_{1} \right) \right\rbrack
```

with the asymptote read off the two anchors on the same powered scale,
$`L_{\infty}^{\rho} = L_{1}^{\rho} + \left( L_{2}^{\rho} - L_{1}^{\rho} \right)/\left( 1 - e^{-K(A_{2} - A_{1})} \right)`$.
Setting $`\rho = 1`$ recovers the von Bertalanffy curve exactly, so the
Richards form nests it and the linear phase below $`A_{1}`$ is
unchanged.

##### Spread of length at age

Fish of one age are spread around the mean length, and that spread grows
with size. Rather than a variance per age, two parameters carry it:
$`CV_{1}`$ and $`CV_{2}`$, the values at the two reference ages, with
linear interpolation between them. Which quantity the interpolation runs
on, and whether the two parameters are read as coefficients of variation
or as standard deviations, are separate choices (`growth_cv_type` and
`growth_sd_type`).

The interpolation. The branches are on real age $`x`$: below $`A_{1}`$
the value is held at $`CV_{1}`$, from $`A_{2}`$ on it is held at
$`CV_{2}`$, and between them it is interpolated. Under
`growth_cv_type = "len"` (the default) the interpolant is mean length,
so a fish whose mean length sits halfway from $`L_{1}`$ to $`L_{2}`$
takes a value halfway from $`CV_{1}`$ to $`CV_{2}`$; under `"age"` it is
age itself:

``` math
c(x) = \begin{cases}
CV_{1} & x < A_{1} \\[1ex]
CV_{1} + \left( CV_{2} - CV_{1} \right)\dfrac{\bar{L}(x) - L_{1}}{L_{2} - L_{1}} & A_{1} \leq x < A_{2}^{cv}, \quad \texttt{cv\_type = "len"} \\[2ex]
CV_{1} + \left( CV_{2} - CV_{1} \right)\dfrac{x - A_{1}}{A_{2}^{cv} - A_{1}} & A_{1} \leq x < A_{2}^{cv}, \quad \texttt{cv\_type = "age"} \\[2ex]
CV_{2} & x \geq A_{2}^{cv}
\end{cases}
```

where $`A_{2}^{cv} = A_{2}`$, or the oldest age being evaluated when
`growth_A2 = "Linf"` and $`L_{2}`$ is read as the asymptote. Note that
the two choices are independent: under `"len"` the branch is still taken
on age, and only the value inside the middle branch is a function of
$`\bar{L}(x)`$. Because $`\bar{L}`$ is itself increasing in $`x`$, the
two give similar shapes; they differ where growth is fast, since
interpolating on length compresses the change into the young ages.

From $`c(x)`$ to the standard deviation. Under `growth_sd_type = "cv"`
(the default) $`c(x)`$ is a coefficient of variation and is multiplied
by mean length, so the spread scales with size; under `"sd"` the two
parameters are standard deviations already and $`c(x)`$ is the standard
deviation, interpolated by the same rule:

``` math
\sigma_{L}(x) = \begin{cases}
c(x)\,\bar{L}(x) & \texttt{sd\_type = "cv"} \\[1ex]
c(x) & \texttt{sd\_type = "sd"}
\end{cases}
```

The two differ in what is held constant across ages. Under `"cv"` a
constant $`CV`$ means the spread grows in proportion to length, so a 20
cm fish and a 60 cm fish have standard deviations in the ratio 1:3;
under `"sd"` a constant value means the same absolute spread at every
age. A $`CV`$ that declines with size ($`CV_{2} < CV_{1}`$) is the usual
estimate, since the lengths of young fish vary with birth date and early
growth while old fish have converged on the asymptote.

##### The plus group

The accumulator age $`a_{+}`$ holds every fish of that age and older, so
its mean length is not the curve at $`a_{+}`$: most of its fish are
older than $`a_{+}`$ and larger. Under `growth_plus_group = "mixture"`
the plus group is treated as a mixture of the ages it contains, with two
assumptions. Its age composition is taken to decline geometrically, the
share of fish $`k`$ years past the accumulator age being proportional to
$`e^{-0.2k}`$, which is the survivorship under a total mortality of
$`0.2`$ per year; the $`0.2`$ is a fixed assumption about how quickly
numbers decline with age, not an estimated rate, and it sets how much
weight the older, larger fish carry. And over those years mean length is
taken to rise linearly from $`\bar{L}(a_{+})`$ to $`L_{\infty}`$ across
a second lifetime of $`a_{+}`$ years, beyond which the curve is at its
asymptote. The plus group’s mean length is the survivorship-weighted
mean of those lengths:

``` math
\bar{L}_{+} = \dfrac{\sum_{k = 0}^{a_{+}}e^{- 0.2k}\left\lbrack \bar{L}(a_{+}) + \dfrac{k}{a_{+}}\left( L_{\infty} - \bar{L}(a_{+}) \right) \right\rbrack}{\sum_{k = 0}^{a_{+}}e^{- 0.2k}}
```

which lies between $`\bar{L}(a_{+})`$ and $`L_{\infty}`$, closer to
$`\bar{L}(a_{+})`$ when $`a_{+}`$ is old enough that little growth
remains. Under `growth_plus_group = "curve"` the plus group is simply
the curve at $`a_{+}`$. Within a year the plus group grows from
$`\bar{L}_{+}`$ by the von Bertalanffy increment over the elapsed time
rather than being re-read from the curve, so its mean length moves
through the season the way every other age’s does.

##### Growth that changes over time

Growth can move over the series in two ways, which answer different
questions and can be used together.

Varying the parameters. Any growth parameter can carry a deviation
series (`growth_tv_model` names one structure per parameter, `"iid"` or
`"rw"`), so the parameter in year $`y`$ is

``` math
P_{y} = \begin{cases}
P\exp(\delta_{y}) & \texttt{growth\_tv\_link = "log"} \\[1ex]
P^{lo} + \left( P^{hi} - P^{lo} \right)\,\mathrm{logit}^{-1}\left\lbrack \mathrm{logit}\!\left( \dfrac{P - P^{lo}}{P^{hi} - P^{lo}} \right) + \delta_{y} \right\rbrack & \texttt{"logit"}
\end{cases}
```

The log link suits a positive parameter. The logit link keeps a
parameter strictly inside bounds $`\left( P^{lo}, P^{hi} \right)`$
however large the deviation, which every growth parameter is. The
deviations are penalized by the same process error used for selectivity
(see Growth Deviations under the penalties), and the realized parameters
are reported year by year as `growth_pars_y`.

How the deviated parameters reach size at age is a second choice
(`growth_tv_type`). Under `"curve"` each year’s sizes are simply read
off that year’s own curve: a fish of age $`a`$ in year $`y`$ has the
length year $`y`$’s parameters give age $`a`$, with no memory of the
years it actually lived through. Under `"cohort"` size at age is carried
forward instead, so each cohort keeps the history of the conditions it
experienced. Writing $`g_{y}`$ for the growth increment over one year
under year $`y`$’s parameters,

``` math
g_{y}(L) = \left\lbrack L_{\infty,y}^{\rho_{y}} + \left( L^{\rho_{y}} - L_{\infty,y}^{\rho_{y}} \right)e^{- K_{y}} \right\rbrack^{1/\rho_{y}}
```

the start-of-year mean length advances as

``` math
\bar{L}_{y+1,a} = \begin{cases}
L_{0} + \dfrac{a}{A_{1}}\left( L_{1,\,y+1-a} - L_{0} \right) & a < A_{1} \quad \text{(the cohort's own birth-year } L_{1}) \\[2ex]
\bar{L}_{y+1}(a) & a = a^{\ast} \quad \text{(the first age past } A_{1}, \text{ read off the year's curve)} \\[2ex]
g_{y}\left( \bar{L}_{y,a-1} \right) & a^{\ast} < a < a_{+} \\[2ex]
\dfrac{\left( N_{y,a_{+}-1} + 0.01 \right)g_{y}\left( \bar{L}_{y,a_{+}-1} \right) + \left( N_{y,a_{+}} + 0.01 \right)g_{y}\left( \bar{L}_{y,a_{+}} \right)}{N_{y,a_{+}-1} + N_{y,a_{+}} + 0.02} & a = a_{+}
\end{cases}
```

Three consequences follow. Ages still in the linear phase take the
length at $`A_{1}`$ their birth year’s parameters gave them, so a cohort
born in a poor year carries that start forward. The plus group is the
only place growth depends on abundance: the cohort just entering it and
the fish already there are blended by their numbers, which is why growth
under this option is evaluated inside the population dynamics year loop
rather than before it. And the coefficient of variation at age is held
at the first year’s, $`c(x)`$ being evaluated once from the first year’s
curve and parameters and then held while the mean moves, so the
deviations change mean size without also changing the spread.

Semi-parametric growth. The second way is a surface of deviations on
mean length at age itself, indexed by year and age (`growth_semipar`),

``` math
\bar{L}_{y,a}^{\ast} = \bar{L}_{y,a}\exp\left( \varepsilon_{y,a} \right)
```

applied after the curve and after any cohort propagation. The parametric
curve stays the parametric part and $`\varepsilon`$ holds departures
from it: a curve with a handful of parameters cannot fit a year in which
only the four-year-olds were small, and an unconstrained transition per
year is not identified. The spread follows the deviated mean through
$`\sigma_{L}`$, so under `growth_cv_type = "len"` a deviation that
lengthens a fish also moves it along the $`CV`$ ramp, while under
`"age"` the spread at age is untouched. Available structures are
`"iid"`, `"rw"`, `"2dar1"` and `"3dmarg"`/`"3dcond"` – the same process
errors the semi-parametric selectivity forms use, described under the
penalties.

##### The key

The key is that distribution integrated over each length bin. Length at
age is normal about $`\bar{L}(a)`$ with standard deviation
$`\sigma_{L}(a)`$, so with $`\ell_{l}`$ the lower edge of bin $`l`$ each
edge is first standardized,

``` math
z_{l,a} = \dfrac{\ell_{l} - \bar{L}(a)}{\sigma_{L}(a)}
```

which is where $`\sigma_{L}`$ enters: it sets how many bins the age’s
mass is spread over, a large $`\sigma_{L}`$ flattening the column and a
small one concentrating it near $`\bar{L}(a)`$. Each bin then takes the
probability between its edges, $`\Phi`$ being the standard normal
distribution function, and the two tails are accumulated into the end
bins so no mass is lost:

``` math
A_{l,a} = \begin{cases}
\Phi(z_{2,a}) & l = 1 \\[1ex]
\Phi(z_{l+1,a}) - \Phi(z_{l,a}) & 1 < l < n_{l} \\[1ex]
1 - \Phi(z_{n_{l},a}) & l = n_{l}
\end{cases}
```

The first bin absorbs everything below $`\ell_{2}`$ and the last
everything above $`\ell_{n_{l}}`$, so the column sums to one,
$`\sum_{l}A_{l,a} = 1`$, for any $`\bar{L}(a)`$ and $`\sigma_{L}(a)`$,
including means that fall outside the binned range entirely. That is
what makes each column a proper conditional distribution $`P(l \mid a)`$
and lets the composition likelihoods normalize.

Under `growth_dist = "lognormal"` the same differences are taken on the
log scale, about a median-corrected mean so that $`E[L] = \bar{L}(a)`$:

``` math
z_{l,a} = \dfrac{\log \ell_{l} - \left( \log \bar{L}(a) - \sigma_{L}^{2}(a)/2 \right)}{\sigma_{L}(a)}
```

Here $`\sigma_{L}(a)`$ is a standard deviation on the log scale, so
`growth_sd_type = "sd"` is the pairing that keeps it dimensionally
consistent; combining `"lognormal"` with `"cv"` feeds a length-scale
spread into a log-scale transform.

##### Weight at age

When `waa_model = "wt_len"`, weight at age is the key applied to weight
at the bin midpoints $`\tilde{\ell}_{l}`$ through the weight-length
relationship $`W = \alpha L^{\beta}`$,

``` math
W_{p,r,y,\tau,a,s} = \sum_{l}A_{p,r,y,\tau,l,a,s}\,\alpha_{p,r,s}\,\tilde{\ell}_{l}^{\beta_{p,r,s}}
```

so the weight at age carries the spread of length at age rather than
being the weight of the mean length.

##### Timing within the year

Every fleet has its own key and weight, read at the point in the season
that fleet’s observations are taken: each fishery fleet’s at `t_fish`,
each survey’s at `t_srv`, and the spawning weight at the spawning time.
A season starts at the cumulative duration of the seasons before it, and
a point inside it is that start plus the fraction elapsed times the
season’s duration. With growth constant over years the curve read at the
real age gives the same mean length whatever the seasons’ durations, so
seasons add no approximation to growth.

#### Fishery Observation Model

The fishery observation model describes the expected retained
catch-at-age, retained catch-at-length, discarded catch-at-age,
discarded catch-at-length, catch and discard (in units of biomass or
abundance), and fishery indices.

Expected retained catch-at-age ($`C_{p,r,y,\tau,a,s,f}^{a}`$) for a
given fishery fleet is calculated using Baranov’s catch equation applied
to the retained fishing mortality:

``` math
\begin{matrix}
C_{p,r,y,\tau,a,s,f}^{a} = \dfrac{\text{retFmort}_{p,r,y,\tau,a,s,f}}{Z_{p,r,y,\tau,a,s}}N_{p,r,y,\tau,a,s}\left\lbrack 1 - \exp\left( - Z_{p,r,y,\tau,a,s} \right) \right\rbrack \\
\end{matrix}
```

Expected dead discarded catch-at-age ($`D_{p,r,y,\tau,a,s,f}^{a}`$) is
similarly:

``` math
\begin{matrix}
D_{p,r,y,\tau,a,s,f}^{a} = \dfrac{\text{discFmort}_{p,r,y,\tau,a,s,f}}{Z_{p,r,y,\tau,a,s}}N_{p,r,y,\tau,a,s}\left\lbrack 1 - \exp\left( - Z_{p,r,y,\tau,a,s} \right) \right\rbrack \\
\end{matrix}
```

These are the forms used under $`\text{move\_timing}`$ 0 and 1, where
individuals occupy a single region for the whole season and catch is
taken where they actually are: at the post-movement (destination)
locations under $`\text{move\_timing} = 0`$, since movement happens at
the start of the season, and at the pre-movement (origin) locations
under $`\text{move\_timing} = 1`$, since it happens at the end. Under
$`\text{move\_timing} = 2`$ fish redistribute among regions while they
are being caught, the region-local Baranov equation is no longer valid,
and catch and discards are computed from the season-integrated abundance
$`N^{\int}`$ instead (see Catch Under Continuous Movement above).

To track length-based dynamics, retained catch-at-length
($`C_{p,r,y,\tau,l,s,f}^{l}`$) and discarded catch-at-length
($`D_{p,r,y,\tau,l,s,f}^{l}`$) are derived using:

``` math
\begin{matrix}
C_{p,r,y,\tau,l,s,f}^{l} = \mathbf{A}_{p,r,y,\tau,s}^{l}{\mathbf{C}^{\mathbf{a}}}_{p,r,y,\tau,s,f} \\
D_{p,r,y,\tau,l,s,f}^{l} = \mathbf{A}_{p,r,y,\tau,s}^{l}{\mathbf{D}^{\mathbf{a}}}_{p,r,y,\tau,s,f} \\
\end{matrix}
```

where $`\mathbf{A}_{p,r,y,\tau,s}^{l}`$ is the size-age transition
matrix, supplied as data or built by the growth module. When conditional
age-at-length data are fit (`do_caal = 1`), the joint retained and
discarded catch at length and age are also formed, each age column of
the transition scaled by the catch at that age:

``` math
\begin{matrix}
C_{p,r,y,\tau,l,a,s,f}^{la} = A_{p,r,y,\tau,l,a,s}\,C_{p,r,y,\tau,a,s,f}^{a} \\
D_{p,r,y,\tau,l,a,s,f}^{la} = A_{p,r,y,\tau,l,a,s}\,D_{p,r,y,\tau,a,s,f}^{a} \\
\end{matrix}
```

so that summing over lengths returns the catch at age and summing over
ages the catch at length.

###### Selecting at length rather than at age

The catch at length above is the catch at age spread over the key, which
is what `FishLenComps_sel = "age"` (the default) and
`SrvLenComps_sel = "age"` give. Selectivity has already been applied at
age by then, as its average over the lengths the age covers,
$`s_{a} = \sum_{l}A_{l,a}\,s_{l}`$. Every fish of an age is therefore
equally catchable, and the length composition within an age is the key’s
own column.

Under `"length"` the order is reversed. The numbers at age are spread
over the key first and selected length by length,

``` math
C_{p,r,y,\tau,l,s,f}^{l} = s_{r,y,l,s,f}\sum_{a}^{a_{+}}A_{p,r,y,\tau,l,a,s}\,N_{p,r,y,\tau,a,s}\left( 1 - e^{- Z_{p,r,y,\tau,a,s}} \right)\frac{F_{r,y,\tau,f}}{Z_{p,r,y,\tau,a,s}}
```

and the survey’s index compositions the same way. The two differ by the
covariance of length and selection within an age. Selecting at age
replaces $`s_{l}`$ by its mean over the age’s length range, which is
exact only where $`s_{l}`$ is flat over that range. Where the curve is
steep, the fish taken from an age are longer or shorter than that age’s
average fish, and one number per age cannot represent that. The
difference is real and not a rounding, so a model whose length
compositions inform a length-based selectivity should select at length.

Expected retained catch ($`\text{Catch}_{r,y,\tau,f}`$) is computed by
summing over populations and then either as abundance:

``` math
\begin{matrix}
\text{Catch}_{r,y,\tau,f} = \sum_{p}^{n_p}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}C_{p,r,y,\tau,a,s,f}^{a}} \\
\end{matrix}
```

or as biomass:

``` math
\begin{matrix}
\text{Catch}_{r,y,\tau,f} = \sum_{p}^{n_p}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}C_{p,r,y,\tau,a,s,f}^{a}}W_{p,r,y,\tau,a,s,f}^{fish} \\
\end{matrix}
```

The fishery weight at age $`W^{fish}`$ is normally the population’s mean
weight at age at the fleet’s timing, which is the weight of an average
fish of that age. The catch is not made of average fish. A length-based
gear takes the long ones from an age more often than the short ones, and
the long ones weigh more. With `fish_waa_selected = 1` the fleet’s catch
biomass instead uses the mean weight of the fish it takes at each age,
the weight averaged over the key re-weighted by selectivity,

``` math
W_{p,r,y,\tau,a,s,f}^{fish} = \dfrac{\sum_{l}A_{p,r,y,\tau,l,a,s}\,s_{r,y,l,s,f}\,w_{l}}{\sum_{l}A_{p,r,y,\tau,l,a,s}\,s_{r,y,l,s,f}}
```

with $`w_{l}`$ the weight at the bin midpoint from `wt_len_pars`. Flat
selectivity returns the population mean exactly, and a knife edge
returns the weight of the one bin it keeps, so the option changes
nothing unless the gear selects within an age. Where it does, the
difference is largest at the youngest selected ages, whose lengths fall
on the ascending limb. `srv_waa_selected` does the same for a survey
index in weight.

Population-specific predicted retained catch
($`\text{Catch}_{p,r,y,\tau,f}`$) retains the population index and is
not summed across $`p`$:

``` math
\begin{matrix}
\text{Catch}_{p,r,y,\tau,f} = \sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}C_{p,r,y,\tau,a,s,f}^{a}}W_{p,r,y,\tau,a,s,f}^{fish} \\
\end{matrix}
```

Expected total discards ($`\text{Discard}_{r,y,\tau,f}`$) are computed
from the dead discarded catch-at-age scaled back by the discard
mortality rate to yield total discarded individuals (dead and released
alive):

``` math
\begin{matrix}
\text{Discard}_{r,y,\tau,f} = \sum_{p}^{n_p}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}\frac{D_{p,r,y,\tau,a,s,f}^{a}}{\delta_{r,y,\tau,f}}} \quad \text{(abundance)} \\
\end{matrix}
```

``` math
\begin{matrix}
\text{Discard}_{r,y,\tau,f} = \sum_{p}^{n_p}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}\frac{D_{p,r,y,\tau,a,s,f}^{a}}{\delta_{r,y,\tau,f}}}W_{p,r,y,\tau,a,s,f}^{fish} \quad \text{(biomass)} \\
\end{matrix}
```

Discards can also be expressed as a fraction of total catch (abundance
or biomass). Population-specific discards follow the same structure
without summing across $`p`$.

Similarly, expected fishery indices ($`\text{FshIdx}_{p,r,y,\tau,f}`$)
can be computed as either abundance-based or biomass-based, using the
product of total fishery selectivity and retention selectivity:

``` math
\begin{matrix}
\text{FshIdx}_{p,r,y,\tau,f} = q_{r,y,f}^{\text{Fsh}}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}N_{p,r,y,\tau,a,s}}\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}}\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}} \\
\end{matrix}
```

where $`q_{r,y,f}^{\text{Fsh}}`$ is the catchability coefficient for a
given fishery fleet. By default the fishery index is computed directly
from start-of-season numbers-at-age without survival discounting, on the
rationale that fishing is spread across the whole season rather than
occurring at a point within it. However, an optional within-season
timing $`t_{r,\tau,f}^{fish}`$ (`t_fish`, a fraction of the season) can
be specified per region, season, and fleet, mirroring the survey
convention: when supplied, the numbers entering the index are first
decayed by

``` math
N_{p,r,y,\tau,a,s}\exp\left( - Z_{p,r,y,\tau,a,s} \cdot t_{r,\tau,f}^{fish} \right)
```

with $`t^{fish} = 0`$ reproducing the start-of-season default. This is
appropriate for fishery CPUE series that reflect a survey-like snapshot
(e.g., an index standardized to a particular part of the season).

The set of ages entering the index sum can also be restricted per fleet
via `fish_idx_ages`, replacing $`\sum_a`$ with a sum over the named ages
only. The restriction applies to the index total alone; the fleet’s
selectivity, catch, and compositions are unaffected.

For the same reason, the fishery index follows catch rather than the
survey when movement is continuous. Under $`\text{move\_timing} = 2`$
the fleet encounters fish as they redistribute, so
$`N_{p,r,y,\tau,a,s}`$ above is replaced by the season-integrated
abundance $`N^{\int}_{p,r,y,\tau,a,s}`$. The integral carries units of
abundance $`\times`$ time, and the resulting constant is absorbed by the
estimated $`q_{r,y,f}^{\text{Fsh}}`$. Under $`\text{move\_timing}`$ 0
and 1 the index uses $`N_{p,r,y,\tau,a,s}`$ as written, at the
post-movement and pre-movement locations respectively. Biomass-based
fishery indices are computed as:

``` math
\begin{matrix}
\text{FshIdx}_{p,r,y,\tau,f} = q_{r,y,f}^{\text{Fsh}}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}N_{p,r,y,\tau,a,s}}\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}}\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}}W_{p,r,y,\tau,a,s,f}^{fish} \\
\end{matrix}
```

The observed region-aggregated fishery index is compared to the sum of
predicted indices across populations:
$`\sum_p \text{FshIdx}_{p,r,y,\tau,f}`$. Population-specific fishery
indices are compared directly to $`\text{FshIdx}_{p,r,y,\tau,f}`$
without summation.

#### Survey Observation Model

Likewise, the survey observation model describes the expected survey
catch-at-age, survey catch-at-length, and survey indices. Expected
survey catch-at-age ($`I_{p,r,y,\tau,a,s,sf}^{a}`$) is calculated as
follows:

``` math
\begin{matrix}
I_{p,r,y,\tau,a,s,sf}^{a} = N_{p,r,y,\tau,a,s}\exp\left( - Z_{p,r,y,\tau,a,s} \cdot t_{r,\tau,sf}^{srv} \right)\text{Sel}_{p,r,y,\tau,a,s,sf}^{\text{Srv}} \\
\end{matrix}
```

where subscript $`sf`$ denotes a given survey fleet,
$`t_{r,\tau,sf}^{srv}`$ is the survey timing as a fraction of the
season, and $`\text{Sel}_{p,r,y,\tau,a,s,sf}^{\text{Srv}}`$ is the
survey selectivity-at-age pattern.

This elementwise discount is exact under $`\text{move\_timing}`$ 0 and
1, because movement for the season has already been resolved and the
fish are stationary within it (at their post-movement locations under
$`\text{move\_timing} = 0`$ and their pre-movement locations under
$`\text{move\_timing} = 1`$). It is not exact under
$`\text{move\_timing} = 2`$, where it would hold fish in place while
they are in fact diffusing. There, the survey observes the population
propagated a fraction $`t^{srv}`$ into the season under the combined
generator $`\mathbf{\Lambda}_{p,y,\tau,a,s}`$ of the Population
Projection section:

``` math
\begin{matrix}
I_{p,r,y,\tau,a,s,sf}^{a} = \left\lbrack \exp\left( \mathbf{\Lambda}_{p,y,\tau,a,s}\,t_{r,\tau,sf}^{srv} \right)\mathbf{N}_{p,y,\tau,a,s} \right\rbrack_{r}\text{Sel}_{p,r,y,\tau,a,s,sf}^{\text{Srv}} \\
\end{matrix}
```

Note that a survey is a snapshot at an instant within the season rather
than an accumulation over it, so this is the partial propagation
$`\exp(\mathbf{\Lambda}\,t^{srv})`$ and not the season integral used for
catch. Selectivity applies at the destination region, i.e. after
propagation.

Survey timing may differ by region, which a single propagation operator
cannot represent: fish observed in region $`r`$ arrived from regions
whose elapsed times differ. The convention adopted is that the survey in
region $`r`$ observes the population propagated to that region’s survey
time, which is what the subscript $`r`$ outside the bracket denotes.
When $`t^{srv}_{r,\tau,sf}`$ is constant across regions, the usual case,
this collapses to a single propagation and one matrix exponential.

Expected survey catch-at-length ($`I_{p,r,y,\tau,l,s,sf}^{l}`$) is given
by:

``` math
\begin{matrix}
I_{p,r,y,\tau,l,s,sf}^{l} = \mathbf{A}_{p,r,y,\tau,s}^{l}{\mathbf{I}^{\mathbf{a}}}_{p,r,y,\tau,s,sf} \\
\end{matrix}
```

and, when conditional age-at-length data are fit, the joint survey index
at length and age
$`I_{p,r,y,\tau,l,a,s,sf}^{la} = A_{p,r,y,\tau,l,a,s}\,I_{p,r,y,\tau,a,s,sf}^{a}`$.

Survey indices ($`\text{SrvIdx}_{p,r,y,\tau,sf}`$) can be computed as
either abundance-based or biomass-based. Abundance-based survey indices
are calculated as:

``` math
\begin{matrix}
\text{SrvIdx}_{p,r,y,\tau,sf} = q_{r,y,sf}^{\text{Srv}}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}I_{p,r,y,\tau,a,s,sf}^{a}} \\
\end{matrix}
```

while biomass-based indices are computed as:

``` math
\begin{matrix}
\text{SrvIdx}_{p,r,y,\tau,sf} = q_{r,y,sf}^{\text{Srv}}\sum_{a}^{a_{+}}{\sum_{s}^{n_{s}}I_{p,r,y,\tau,a,s,sf}^{a}}W_{p,r,y,\tau,a,s,sf}^{srv} \\
\end{matrix}
```

Here, $`W_{p,r,y,\tau,a,s,sf}^{srv}`$ is the weight-at-age for a given
survey, $`q_{r,y,sf}^{\text{Srv}}`$ represents the survey catchability
coefficient. The observed region-aggregated survey index is compared to
the sum of predicted indices across populations:
$`\sum_p \text{SrvIdx}_{p,r,y,\tau,sf}`$. Population-specific survey
indices are compared directly to $`\text{SrvIdx}_{p,r,y,\tau,sf}`$
without summation.

A survey fleet can instead observe year class strength directly, with
`srv_idx_type = "recdev"`. Such a fleet reads no part of the population:

``` math
\begin{matrix}
\text{SrvIdx}_{p,r,y,\tau,sf} = q_{r,y,sf}^{\text{Srv}}\left( \varepsilon^{R}_{p,r,y} - \mu^{R}_{y} \right) \\
\end{matrix}
```

where $`\varepsilon^{R}_{p,r,y}`$ is the recruitment deviation and
$`\mu^{R}_{y}`$ the center its penalty asserts (see Recruitment). The
deviation relative to that center is the anomaly, how strong the year
class was against what the model expected, which is what a pre-recruit
survey or an environmental index measures. Under a bias ramp
$`\mu^{R}_{y} = -b_{y}\sigma_{R}^{2}/2`$, so the anomaly and the
deviation are not the same quantity. Deviations are signed, so such a
fleet requires a normal index likelihood, and its selectivity, survey
timing and weight-at-age are never read.

As with the fishery index, the ages entering the survey index sum can be
restricted per fleet via `srv_idx_ages`. Restricting a fleet to a single
age turns it into an index of that age alone (e.g., an age-1 acoustic
index of recruitment strength), while the fleet’s compositions continue
to use the full age range because the restriction applies to the index
sum rather than to selectivity.

Survey catchability is by default an estimated parameter,
$`q_{r,y,sf}^{\text{Srv}} = \exp\left(\ln q_{r,\hat{b}(y),sf}\right)`$
with time-block structure $`\hat{b}(y)`$. Alternatively, catchability
can be concentrated out of the likelihood analytically (`srv_q_type`),
treating it as a pure scaling nuisance parameter. Two analytic solutions
are available, computed per region and fleet over only the years with
observations. The arithmetic solution (`"arith"`) is the ratio of mean
observed to mean predicted unscaled index:

``` math
\hat{q}_{r,sf} = \frac{\overline{\text{ObsSrvIdx}_{r,\cdot,sf}}}{\overline{\text{SrvIdx}_{r,\cdot,sf}^{(q=1)}}}
```

and the geometric solution (`"geo"`) is the exponentiated mean
log-ratio:

``` math
\hat{q}_{r,sf} = \exp\left( \overline{\log\left(\text{ObsSrvIdx}_{r,\cdot,sf}\right) - \log\left(\text{SrvIdx}_{r,\cdot,sf}^{(q=1)}\right)} \right)
```

where $`\text{SrvIdx}^{(q=1)}`$ denotes the predicted index evaluated at
$`q = 1`$. The geometric form is the exact maximum likelihood solution
for a lognormal index likelihood with a shared standard error, and is
the usual companion to a lognormal index; the arithmetic form matches
the convention some existing assessments use. Analytic fleets fix their
$`\ln q^{\text{Srv}}`$ parameters automatically, ignore any block
structure, and cannot carry catchability covariates or priors.

For estimated survey catchability, environmental linkage can be
specified:

``` math
q_{r,y,sf}^{\text{Srv}} = q_{r,sf}^{\text{Srv}}\exp\left( \mathbf{x}^{T}\mathbf{\beta +}\sum_{m}^{}{\iota_{m}p_{m}\left( z_{r,y,sf} \right)} \right)
```

where $`q_{r,sf}^{\text{Srv}}`$ is the base survey catchability (i.e.,
intercept), $`\mathbf{x}`$ is a matrix of covariates, $`\mathbf{\beta}`$
is a vector of regression coefficients, $`\iota_{m}p_{m}`$ are
orthogonal polynomial coefficients along with its basis functions, and
$`z_{r,y,sf}`$ are the covariates for which a polynomial term is
assumed.

#### Tagging Observation Model

The tagging observation model tracks tag cohorts
($`T_{p,r,y,\tau,a,s}^{k}`$) by the combination of release region,
release year, and release season ($`k`$) and follows a Brownie tag
attrition framework. Tag cohorts are tracked for a pre-defined maximum
duration (maximum tag liberty; $`n_{L}`$), after which calculations for
the tag cohort are no longer computed. Tag dynamics incorporate both
population ($`p`$) and season ($`\tau`$) dimensions, and tag reporting
rates are fleet-specific ($`\beta_{r,y,f}`$). In general, the process
dynamics for the tagged cohort mimic those specified for the overall
population. Immediately following release, tag cohorts are decremented
by an initial tag-induced mortality rate:

``` math
T_{p,r,y,\tau,a,s}^{k} = T_{p,r,y,\tau,a,s}^{k}\exp( - \eta^{\text{mort}})
```

where $`\eta^{\text{mort}}`$ is the initial tag-induced mortality rate.

Tag cohorts are then advanced across each season by the same seasonal
transition operator as the population, but evaluated with the
tag-specific total mortality $`\mathbf{Z}^{\text{Tag}}`$ and over the
fraction of the season the cohort is actually at liberty for,
$`\Delta\tau\,\omega^{k}_{\tau}`$ (with $`\omega^{k}_{\tau}`$ defined
below):

``` math
\mathbf{\Phi}_{p,y,\tau,a,s}^{\text{Tag}} = \begin{cases}
\mathbf{M}_{p,y,\tau,a,s}\,\text{diag}\left( \exp\left( -\mathbf{Z}^{\text{Tag}}_{p,y,\tau,a,s} \right) \right) & \text{move\_timing} = 0 \\[4pt]
\text{diag}\left( \exp\left( -\mathbf{Z}^{\text{Tag}}_{p,y,\tau,a,s} \right) \right)\mathbf{M}_{p,y,\tau,a,s} & \text{move\_timing} = 1 \\[4pt]
\left[ \exp\left( \dot{\mathbf{Q}}_{p,y,\tau,a,s}^{T}\,\Delta\tau\,\omega^{k}_{\tau} - \text{diag}\left( \mathbf{Z}^{\text{Tag}}_{p,y,\tau,a,s} \right) \right) \right]^{T} & \text{move\_timing} = 2 \\
\end{cases}
```

For seasons within a year ($`\tau < n_\tau`$):

``` math
\left( \mathbf{T}_{p,y,\tau + 1,a,s}^{k} \right)^{T} = \left( \mathbf{T}_{p,y,\tau,a,s}^{k} \right)^{T}\mathbf{\Phi}_{p,y,\tau,a,s}^{\text{Tag}}
```

and at the end of the final season ($`\tau = n_\tau`$) individuals
advance in age, with accumulation of individuals in the plus group,
exactly as for the population:

``` math
\left( \mathbf{T}_{p,y + 1,1,a + 1,s}^{k} \right)^{T} = \left( \mathbf{T}_{p,y,n_\tau,a,s}^{k} \right)^{T}\mathbf{\Phi}_{p,y,n_\tau,a,s}^{\text{Tag}},\quad\text{for }1 \leq a < a_{+}
```

``` math
\left( \mathbf{T}_{p,y + 1,1,a_{+},s}^{k} \right)^{T} = \left( \mathbf{T}_{p,y+1,1,a_{+},s}^{k} \right)^{T} + \left( \mathbf{T}_{p,y,n_\tau,a_{+},s}^{k} \right)^{T}\mathbf{\Phi}_{p,y,n_\tau,a_{+},s}^{\text{Tag}}
```

Under $`\text{move\_timing} = 0`$ these reduce to the historical
sequential form, in which Markovian movement
$`\left( \mathbf{T}_{p,y,\tau,a,s}^{k} \right)^{\text{T}}\mathbf{M}_{p,y,\tau,a,s}`$
is applied first and the cohort then follows an elementwise exponential
mortality model,
$`T_{p,r,y,\tau + 1,a,s}^{k} = T_{p,r,y,\tau,a,s}^{k}\exp\left( - Z_{p,r,y,\tau,a,s}^{\text{Tag}} \right)`$.

Mid-season releases are handled differently by the two families of
timings. Under the sequential timings, a cohort released partway through
its release season ($`t^{\text{tag}} < 1`$) skips the discrete movement
step entirely for that season, $`\mathbf{M}`$ is replaced by the
identity, because a full-season transition matrix cannot represent a
partial interval. Such a cohort stays in its release region for the
remainder of the release season, experiencing only the partial-interval
mortality, and movement resumes normally from the following season.

Continuous movement needs no such exemption, and imposing one would be
inconsistent: the generator is scaled by
$`\Delta\tau\,\omega^{k}_{\tau}`$, the same at-liberty fraction that
already scales the cohort’s mortality, so tags released mid-season
diffuse for exactly the fraction of the season they were at liberty for.
Freezing them in place while still discounting their mortality partially
would have tags dying on a partial season but holding station for a
whole one.

Total mortality for the tagged cohort
($`Z_{p,r,y,\tau,a,s}^{\text{Tag}}`$) is specified as:

``` math
Z_{p,r,y,\tau,a,s}^{\text{Tag}} = \omega^{k}_{\tau} \left( \kappa \cdot \Delta\tau + \text{NatMort}_{p,r,y,a,s} \cdot \Delta\tau + \sum_{f \in \mathcal{F}^{\text{Tag}}} \left[\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}} \cdot \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}} \cdot \text{Fmort}_{r,y,\tau,f} + \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}} \cdot \left(1 - \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}}\right) \cdot \delta_{r,y,\tau,f} \cdot \text{Fmort}_{r,y,\tau,f}\right] \right)
```

where $`\omega^{k}_{\tau}`$ is the fraction of season $`\tau`$ that
cohort $`k`$ is actually at liberty for. It equals
$`t^{\text{tag}}_{k}`$ in the release season of the release year, and
$`1`$ in every subsequent season:

``` math
\omega^{k}_{\tau} = \begin{cases} t^{\text{tag}}_{k} & \text{release season of the release year} \\ 1 & \text{otherwise} \end{cases}
```

Note that $`\omega^{k}_{\tau}`$ multiplies every mortality component,
fishing included, rather than only the total. This matters for the
recapture equation below: $`F/Z^{\text{Tag}}`$ is the fraction of deaths
owing to fishing and must lie in $`[0,1]`$. Scaling $`Z^{\text{Tag}}`$
while leaving $`F`$ at full-season scale would give
$`F/(Z^{\text{Tag}} t^{\text{tag}})`$, which exceeds $`1`$ whenever
$`t^{\text{tag}} < F/Z^{\text{Tag}}`$ and would predict more recaptures
than there are dead tags. Scaling both leaves the ratio at
$`F/Z^{\text{Tag}}`$ and lets the $`1 - \exp(-Z^{\text{Tag}})`$ term
carry the shorter exposure, which is the standard partial-interval form
of Baranov’s equation.

Here $`\kappa`$ is a parameter describing chronic tag loss (i.e., annual
tag shedding) and $`\mathcal{F}^{\text{Tag}}`$ denotes the subset of
fishing fleets that contribute tagging data. The summation over
$`\mathcal{F}^{\text{Tag}}`$ rather than all fleets is intentional:
restricting the tag mortality calculation to fleets with tagging data
prevents non-tagging fleets from unintentionally influencing tag-based
parameter estimates (e.g., selectivity, reporting rates).

Similar to computations for retained catch-at-age, tag recaptures are
calculated using a modified version of Baranov’s catch equation, with
fleet-specific tag reporting rates applied to the retained component:

``` math
\text{Recap}_{p,r,y,\tau,a,s,f}^{k} = \beta_{r,y,f}\dfrac{\text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}} \cdot \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}} \cdot \text{Fmort}_{r,y,\tau,f}}{Z_{p,r,y,\tau,a,s}^{\text{Tag}}}T_{p,r,y,\tau,a,s}^{k}\left\lbrack 1 - \exp\left( - Z_{p,r,y,\tau,a,s}^{\text{Tag}} \right) \right\rbrack
```

As for the fishery, this region-local form applies under
$`\text{move\_timing}`$ 0 and 1. Under $`\text{move\_timing} = 2`$ tags
redistribute among regions while they are being caught, so recaptures
use the season-integrated tag abundance:

``` math
\mathbf{T}_{p,y,\tau,a,s}^{k,\int} = \left\lbrack \int_{0}^{1}\exp\left( \mathbf{\Lambda}^{\text{Tag}}_{p,y,\tau,a,s}\upsilon \right)d\upsilon \right\rbrack\mathbf{T}_{p,y,\tau,a,s}^{k},\qquad
\mathbf{\Lambda}^{\text{Tag}}_{p,y,\tau,a,s} = \dot{\mathbf{Q}}_{p,y,\tau,a,s}^{T}\Delta\tau\,\omega^{k}_{\tau} - \text{diag}\left( \mathbf{Z}^{\text{Tag}}_{p,y,\tau,a,s} \right)
```

``` math
\text{Recap}_{p,r,y,\tau,a,s,f}^{k} = \beta_{r,y,f} \cdot \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Fsh}} \cdot \text{Sel}_{p,r,y,\tau,a,s,f}^{\text{Ret}} \cdot \text{Fmort}_{r,y,\tau,f} \cdot T_{p,r,y,\tau,a,s}^{k,\int}
```

The integral is evaluated by the same block-matrix construction used for
the population, with $`\mathbf{\Lambda}^{\text{Tag}}`$ in place of
$`\mathbf{\Lambda}`$ (see Evaluating the season integral above). Because
both $`\mathbf{Z}^{\text{Tag}}`$ and the fishing mortality in the
numerator carry the same $`\omega^{k}_{\tau}`$ scaling, this integrates
fishing mortality over exactly the at-liberty fraction of the season:
the substitution $`u = \omega^{k}_{\tau}\upsilon`$ turns the expression
into
$`\int_{0}^{\omega^{k}_{\tau}}\mathbf{F}\exp\left( \mathbf{\Lambda}u \right)\mathbf{T}^{k}du`$
on the unscaled full-season generator, which is the continuous-movement
counterpart of the partial-interval Baranov equation described above.

Here $`\beta_{r,y,f}`$ represents a fleet-specific tag reporting rate
parameter that can vary across space, time, and fleet, which is
estimated in logit space such that it is constrained between $`[0,1]`$.
Recaptures are computed only for fleets in $`\mathcal{F}^{\text{Tag}}`$.

### Fishery and Survey Selectivity

In the following descriptions of selectivity, we omit subscripts for
sexes and fleets for brevity, although note that the equations remain
specific to those model partitions. Several approaches are available for
parameterizing fishery and survey selectivity. Selectivity can be
defined as either age- or length-based. Selectivity parameters are
estimated by region ($`r`$), year ($`y`$), sex ($`s`$), and fleet, and
are explicitly invariant across populations ($`p`$) and seasons
($`\tau`$).

For age-based selectivity, an age vector is applied directly with a
chosen functional form, and the resulting selectivity-at-age
($`\text{Sel}_{p,r,y,\tau,a,s,f}`$) is likewise invariant across
populations and seasons. For length-based selectivity, a length vector
is used to compute selectivity-at-length ($`\text{Sel}^l_{r,y,s,f}`$),
which is then converted to selectivity-at-age via a dot product with the
size-age transition matrix:

``` math
\mathbf{Sel}^{a}_{p,r,y,\tau,s,f} = \left( \mathbf{A}_{p,r,y,\tau,s}^{l} \right)^{T} \mathbf{Sel}^{l}_{r,y,s,f}
```

Because the size-age transition matrix $`\mathbf{A}_{p,r,y,\tau,s}^{l}`$
varies across populations and seasons, the derived selectivity-at-age
inherits population and season specificity upon conversion, even though
the underlying selectivity-at-length parameters remain shared across
these dimensions. Given the age-based nature of the model,
selectivity-at-age is utilized for all subsequent calculations. In the
following we use the subscript $`b`$ to denote a generalized bin number.

Two forms of logistic selectivity can be specified. The first form is
defined as:

``` math
\begin{matrix}
{Sel}_{b} = \frac{1}{1 + \exp\left\lbrack - k\left( b - b^{50} \right) \right\rbrack} \\
\end{matrix}
```

where $`k`$ determines the slope/steepness of the logistic curve and
$`b^{50}`$ is the bin-at-50% selection. Parameters are supplied in the
order $`\left( b^{50},k \right)`$, and every selectivity parameter is
estimated on the log scale, so a starting value is given as $`\log`$ of
the quantity above. The second form can be expressed as:

``` math
\begin{matrix}
{Sel}_{b} = \frac{1}{1 + 19^{\left( \frac{b^{50} - b}{b^{95}} \right)}} \\
\end{matrix}
```

Here, $`b^{50}`$ is also the bin-at-50% selection and $`b^{95}`$ is the
bin-at-95% selection. Beyond the specification of flat-topped
selectivity, `SPoRC` also allows for dome-shaped selectivity. In
particular, a reparametrized gamma function can be specified:

``` math
\begin{matrix}
p = 0.5\left\lbrack \sqrt{\left( b^{\max} \right)^{2} + 4\gamma^{2}} - b^{\max} \right\rbrack \\
\end{matrix}
```

``` math
{Sel}_{b} = \left( \frac{b}{b^{\max}} \right)^{\frac{b^{\max}}{p}}\exp\left( \frac{b^{\max} - b}{p} \right)
```

In this parameterization, $`p`$ is a derived power parameter, $`\gamma`$
is the shape parameter that describes the steepness of the descending
limb, and $`b^{\max}`$ describes the bin-at-maximum selection.
Dome-shaped selectivity can also be expressed as a power function:

``` math
\begin{matrix}
{Sel}_{b} = \frac{1}{b^{\phi}} \\
\end{matrix}
```

with $`\phi`$ being a power parameter that determines the descending
limb of the curve (larger values are steeper). The last dome-shaped
selectivity form that can be specified includes a 6 parameter (denoted
as $`{\widehat{p}}_{1}`$ through $`{\widehat{p}}_{6}`$) double normal
functional form with the following transformations applied:

``` math
p_{1} = \ {\widehat{p}}_{1}
```

``` math
p_{2} = \ p_{1} + w + \frac{0.99 \cdot max(b) - p_{1} - w}{1 + exp( - {\widehat{p}}_{2})}
```

``` math
p_{3} = \exp\left( {\widehat{p}}_{3} \right),\ \ p_{4} = \exp\left( {\widehat{p}}_{4} \right)
```

``` math
p_{5} = \frac{1}{1 + \exp\left( - {\widehat{p}}_{5} \right)},\ \ p_{6} = \frac{1}{1 + \exp\left( - {\widehat{p}}_{6} \right)}
```

$`p_{1}`$ is the bin at which the plateau begins, carried on the bin
scale rather than transformed so that it means a bin, $`p_{2}`$ the bin
at which the plateau ends, $`w`$ the bin width, $`p_{3}`$ and $`p_{4}`$
control the ascending and descending widths, and $`p_{5}`$ and $`p_{6}`$
are the selectivity at the first and last bins. Writing $`A`$ for the
height the two limbs are built up to and the plateau sits at, which is
one unless a sex carries an apical offset (see Sex Offsets below), the
curve is assembled from

``` math
\begin{matrix}
{asc}_{b} = \exp\left( - \frac{\left( b - p_{1} \right)^{2}}{p_{3}} \right),\ \ {asc}_{min} = \exp\left( - \frac{\left( min(b) - p_{1} \right)^{2}}{p_{3}} \right) \\
{ascscaled}_{b} = p_{5} + \left( A - p_{5} \right) \cdot \frac{{asc}_{b} - {asc}_{min}}{1 - {asc}_{min}} \\
{desc}_{b} = \exp\left( - \frac{\left( b - p_{2} \right)^{2}}{p_{4}} \right),\ \ {desc}_{min} = \exp\left( - \frac{\left( max(b) - p_{2} \right)^{2}}{p_{4}} \right) \\
{descscaled}_{b} = A + \left( p_{6} - A \right) \cdot \frac{{desc}_{b} - 1}{{desc}_{min} - 1} \\
join1_{b} = \frac{1}{1 + \exp\left( - 20 \cdot \frac{b - p_{1}}{1 + \left| b - p_{1} \right|} \right)} \\
join2_{b} = \frac{1}{1 + \exp\left( - 20 \cdot \frac{b - p_{2}}{1 + \left| b - p_{2} \right|} \right)} \\
\end{matrix}
```

and are joined together as:

``` math
\begin{matrix}
{Sel}_{b} = {ascscaled}_{b} \cdot \left( 1 - join1_{b} \right) + join1_{b} \cdot \left\lbrack A \cdot \left( 1 - join2_{b} \right) + {descscaled}_{b} \cdot join2_{b} \right\rbrack \\
\end{matrix}
```

Each limb is rescaled by its own value at the bin its endpoint parameter
refers to, which is what makes those parameters mean what they are
named: dividing the ascending limb by $`{asc}_{min}`$ puts $`{Sel}`$ at
exactly $`p_{5}`$ in the first bin, and dividing the descending limb by
$`{desc}_{min}`$ puts it at exactly $`p_{6}`$ in the last. Both anchors
are the model’s own first and last bin, so they follow the bin range
rather than sitting at a fixed value.

The double normal functional form is incredibly flexible and is able to
reduce to both flat-topped and dome-shaped selectivity forms, depending
on the values of the parameters.

Non-parametric selectivity can also be specified, where bin-specific
logit-scale parameters are transformed via the logistic function:

``` math
\begin{matrix}
{Sel}_{b} = \text{logistic}\left( \eta_{b} \right)
\end{matrix}
```

where $`\eta_{b}`$ is a freely estimated logit-scale selectivity
parameter for each bin.

A second non-parametric form (`"nonparlog"`) holds the free parameters
on the log scale instead and standardizes within each year so
selectivity averages to one over a set of bins $`\mathcal{B}`$:

``` math
\begin{matrix}
{Sel}_{y,b} = \frac{\exp\left( \eta_{y,b} \right)}{\frac{1}{\left| \mathcal{B} \right|}\sum_{b' \in \mathcal{B}}^{}{\exp\left( \eta_{y,b'} \right)}}
\end{matrix}
```

The standardization window $`\mathcal{B}`$ is given per fleet
(`fish_sel_norm_bins` and `srv_sel_norm_bins`) and is every bin by
default, i.e. $`\mathcal{B} = \left\{ 1,\ldots,n_{b} \right\}`$ and the
divisor is the mean over the whole curve. It is expressed in whichever
bin domain the fleet’s selectivity is defined on, ages or lengths.
Standardization is applied to every bin of the curve regardless of which
bins $`\mathcal{B}`$ holds; the window sets what the average is taken
over, not which bins are rescaled. Retention selectivity takes the whole
bin range, since the window is supplied for the fishery and survey
streams only.

The two non-parametric forms differ in both respects: the logit form
bounds every raw value below one via the logistic transform and mean
standardizes over years and bins jointly, whereas the log form leaves
the scale free and centers within the year, the convention several
existing assessments use. Under the log form only the differences among
$`\eta`$ within a year are identified; the level of $`\eta`$ is free and
is absorbed by catchability or fishing mortality, which is why it is
typically paired with the selectivity parameter centering penalty
described in the Priors and Penalties section.

Which bins $`\mathcal{B}`$ holds is a statement about the gear rather
than a numerical convenience, because the standardization fixes the
level that catchability is defined against: a fleet whose $`q`$ refers
to only part of the bin range standardizes over that part. Two choices
of $`\mathcal{B}`$ differ only by a constant multiplier on
$`{Sel}_{y,b}`$, and a freely estimated $`q`$ absorbs that multiplier
exactly, leaving every fitted quantity unchanged. Under an informative
prior on $`q`$ that multiplier is not free. The prior resists the
compensating shift in $`\log q`$, and the difference between windows is
carried by the selectivity parameters instead, levered by the derivative
of the index negative log-likelihood with respect to $`\log q`$ at the
fitted values. A window that does not match the definition the
catchability prior was built under therefore shows up as a standing
gradient on $`\eta`$ rather than as an absorbed scalar, and the two
windows give genuinely different fits.

Two additional logistic forms with a freely estimated asymptote
parameter $`\alpha \in (0, 1)`$ are also available. The first uses the
bin-at-50% and slope parameterization:

``` math
\begin{matrix}
{Sel}_{b} = \frac{\alpha}{1 + \exp\left\lbrack - k\left( b - b^{50} \right) \right\rbrack}
\end{matrix}
```

The second uses the bin-at-50% and bin-at-95% parameterization:

``` math
\begin{matrix}
{Sel}_{b} = \frac{\alpha}{1 + 19^{\left( \frac{b^{50} - b}{b^{95}} \right)}}
\end{matrix}
```

where $`\alpha`$ is estimated on the logit scale, allowing the
asymptotic selectivity to be less than 1. These forms are useful for
fleets where full vulnerability is not achieved even at the largest
observed sizes or ages.

##### Selectivity Plateau

Any functional form can additionally carry a plateau (the
`_NSelBins_<n>` suffix on the model string): bins beyond a chosen bin
$`n^{\text{sel}}`$ are held at that bin’s computed value rather than
evaluated through the form,

``` math
{Sel}_{b} = {Sel}_{n^{\text{sel}}} \quad \text{for } b > n^{\text{sel}}
```

applied after the form and its parameter deviations. This is the plateau
convention many existing assessments use, and it is part of the model
rather than a display choice whenever the curve has not saturated by
$`n^{\text{sel}}`$.

##### Sex Offsets

For models with more than one sex, a fleet’s sexes can be linked through
offsets rather than estimated independently or forced identical
(`fish_sel_sex_offset` / `ret_sel_sex_offset` / `srv_sel_sex_offset`,
for total fishery, retention, and survey selectivity respectively).
Under a parameter offset (`"par"`), the stored parameter slots of every
sex beyond the first hold additive offsets on the first sex’s
transformed-scale parameters,

``` math
\theta_{s} = \theta_{1} + \delta_{s}
```

so that for log-scale parameters the second sex’s natural value is the
first sex’s times $`e^{\delta_{s}}`$ (e.g.,
$`k_{M} = k_{F}e^{\delta_{k}}`$,
$`b_{M}^{50} = b_{F}^{50}e^{\delta_{50}}`$). Offsets fixed at zero
reproduce sex-shared parameters. Under a scale offset (`"scale"`), each
sex keeps its own parameters and the realized curve of every sex beyond
the first is multiplied by a constant,

``` math
{Sel}_{b,s} = {Sel}_{b,s}^{\text{form}}\,e^{\gamma_{s}}
```

with $`\gamma_{s}`$ estimated per region, block, and sex; the scaled
curve may exceed one, matching the convention of a log-scale male
selectivity offset applied to the whole curve. Because a constant
multiplier is canceled by mean standardization, the scale offset is
refused for the non-parametric forms and the semi-parametric
time-varying structures.

Under an apical offset (`"apical"`), available for the double normal
alone, the same estimated constant

``` math
A_{s} = e^{\gamma_{s}}
```

is used in a different place. Rather than multiplying the finished
curve, it replaces the one that the double normal’s limbs are built up
to. Writing the two side by side for a sex $`s`$ beyond the first, with
$`{asc}_{b}`$, $`{desc}_{b}`$, $`{asc}_{min}`$, $`{desc}_{min}`$,
$`join1_{b}`$ and $`join2_{b}`$ exactly as defined for the double normal
above, a scale offset evaluates the form with its plateau at one and
then multiplies,

``` math
\begin{matrix}
{ascscaled}_{b} = p_{5} + \left( 1 - p_{5} \right)\dfrac{{asc}_{b} - {asc}_{min}}{1 - {asc}_{min}},\ \ \ 
{descscaled}_{b} = 1 + \left( p_{6} - 1 \right)\dfrac{{desc}_{b} - 1}{{desc}_{min} - 1} \\
\\
{Sel}_{b,s} = \underbrace{e^{\gamma_{s}}}_{\text{outside}} \cdot \left\lbrack {ascscaled}_{b}\left( 1 - join1_{b} \right) + join1_{b}\left( 1 \cdot \left( 1 - join2_{b} \right) + {descscaled}_{b}\,join2_{b} \right) \right\rbrack \\
\end{matrix}
```

while an apical offset puts $`A_{s}`$ everywhere that one appeared,

``` math
\begin{matrix}
{ascscaled}_{b,s} = p_{5} + \left( \underbrace{A_{s}}_{\text{inside}} - p_{5} \right)\dfrac{{asc}_{b} - {asc}_{min}}{1 - {asc}_{min}},\ \ \ 
{descscaled}_{b,s} = \underbrace{A_{s}}_{\text{inside}} + \left( p_{6} - A_{s} \right)\dfrac{{desc}_{b} - 1}{{desc}_{min} - 1} \\
\\
{Sel}_{b,s} = {ascscaled}_{b,s}\left( 1 - join1_{b} \right) + join1_{b}\left( \underbrace{A_{s}}_{\text{inside}} \cdot \left( 1 - join2_{b} \right) + {descscaled}_{b,s}\,join2_{b} \right) \\
\end{matrix}
```

and the curve is never multiplied by anything afterwards. The
consequence is at the two ends. Both arrangements put the plateau at
$`e^{\gamma_{s}}`$, but a scale offset carries $`p_{5}`$ and $`p_{6}`$
down with it while an apical offset leaves them alone:

|  | first bin | plateau | last bin |
|----|----|----|----|
| no offset | $`p_{5}`$ | $`1`$ | $`p_{6}`$ |
| scale offset | $`e^{\gamma_{s}}p_{5}`$ | $`e^{\gamma_{s}}`$ | $`e^{\gamma_{s}}p_{6}`$ |
| apical offset | $`p_{5}`$ | $`A_{s} = e^{\gamma_{s}}`$ | $`p_{6}`$ |

Either can be made to draw any single curve, by rescaling $`p_{5}`$ and
$`p_{6}`$ by $`e^{\gamma_{s}}`$; they differ in how the curve responds
when $`\gamma_{s}`$ moves, which is what a gradient sees.

The apical offset is restricted to the double normal because only there
is it a distinct operation. For a form that peaks at one with no
separate parameters at its ends, scaling and lowering the peak are the
same thing: the ascending logistic scaled by $`A_{s}`$ is the asymptotic
logistic with asymptote $`A_{s}`$ (forms 6 and 7 above), and the gamma
dome is already zero at both ends.

Time-varying deviations are unaffected by any of these options and still
apply per sex to the effective parameters.

##### Bicubic Spline Selectivity

Rather than a fixed functional form, selectivity can instead be
constructed as a smooth two-dimensional surface over bins and years
using a bicubic natural cubic spline. A sparse grid of
$`n_{\hat{b}} \times n_{\hat{y}}`$ freely estimated log-scale node
parameters, $`\eta_{\hat{y},\hat{b}}`$, indexed by bin-node
$`\hat{b} = 1,\ldots,n_{\hat{b}}`$ and year-node
$`\hat{y} = 1,\ldots,n_{\hat{y}}`$, is expanded to the full bin-by-year
surface via two successive natural cubic spline interpolations.

Bin-nodes and evaluation bins are first placed on a common
$`\lbrack 0,1 \rbrack`$ scale, equally spaced by index:

``` math
\tilde{b}_{\hat{b}} = \frac{\hat{b} - 1}{n_{\hat{b}} - 1},\qquad \tilde{b}_{b} = \frac{b - 1}{n_{b} - 1}
```

and a natural cubic spline is fit through the $`n_{\hat{b}}`$ node
positions $`\tilde{b}_{\hat{b}}`$, producing an
$`n_{b} \times n_{\hat{b}}`$ interpolation weight matrix
$`\mathbf{W}^{\text{bin}}`$ such that, for any vector of node values,
$`\mathbf{W}^{\text{bin}}`$ maps those node values onto all $`n_{b}`$
evaluation bins while passing exactly through the node values
themselves. An analogous weight matrix $`\mathbf{W}^{\text{yr}}`$
($`n_{y} \times n_{\hat{y}}`$) is constructed for the year dimension
using year-nodes similarly placed on $`\lbrack0,1\rbrack`$. The full
surface is then obtained by a two-pass tensor-product spline: first
interpolating across bins for every year-node,

``` math
\mathbf{\Xi} = \mathbf{N}\left(\mathbf{W}^{\text{bin}}\right)^{T}
```

where $`\mathbf{N}`$ is the $`n_{\hat{y}} \times n_{\hat{b}}`$ matrix of
node parameters
($`\mathbf{N}_{\hat{y},\hat{b}} = \eta_{\hat{y},\hat{b}}`$) and
$`\mathbf{\Xi}`$ is $`n_{\hat{y}} \times n_{b}`$, and then interpolating
the resulting bin-interpolated year-node curves across years for a given
year $`y`$:

``` math
\log\left(\text{Sel}_{y,b}\right) = \mathbf{W}_{y,}^{\text{yr}}\mathbf{\Xi}_{,b}
```

so that
$`\text{Sel}_{y,b} = \exp\left(\log\left(\text{Sel}_{y,b}\right)\right)`$
for every bin $`b = 1,\ldots,n_{b}`$. Setting $`n_{\hat{y}} = 1`$
collapses the year dimension to a single node (equal weight $`1`$ for
every year), yielding a time-invariant bin-only spline; combining
$`n_{\hat{y}} = 1`$ with discrete time blocks (see Temporal Variation
below) re-fits an independent bin-only spline within each block.

Two optional restrictions can be applied to the range over which the
surface is actually spline-fit, with everything outside that range held
constant (edge-held) rather than continuing the spline.

The first restricts the year dimension: given a user-specified calendar
year $`y^{\text{SelStyr}}`$ within a given block, only years from
$`y^{\text{SelStyr}}`$ through the block’s final year are used to place
year-nodes and evaluate the spline. Years within the block prior to
$`y^{\text{SelStyr}}`$ are assigned the same interpolation weights as
$`y^{\text{SelStyr}}`$ itself (the boundary node),

``` math
\text{Sel}_{y,b} = \text{Sel}_{y^{\text{SelStyr}},b},\qquad y < y^{\text{SelStyr}}
```

i.e., “filled” forward from the first actually-fitted year.

The second restricts the bin dimension: given a user-specified number of
bins $`n_{b}^{\text{fit}} \leq n_{b}`$, bin-nodes and the spline are
only evaluated over bins $`1,\ldots,n_{b}^{\text{fit}}`$; any remaining
bins are held at the last fitted bin’s value,

``` math
\text{Sel}_{y,b} = \text{Sel}_{y,n_{b}^{\text{fit}}},\qquad b > n_{b}^{\text{fit}}
```

This is useful, for example, when the observed age or length range used
to originally fit the surface is narrower than the full number of ages
or lengths represented in the population dynamics.

In addition to the functional forms that can be specified to describe
selectivity processes, several options exist to specify continuous
time-varying processes. In particular, options to specify time-varying
parametric selectivity and time-varying semi-parametric selectivity are
available. To illustrate, if logistic selectivity is specified and
parametric deviations are invoked, the following expression is used:

``` math
{\begin{matrix}
{Sel}_{y,b} = \frac{1}{1 + \exp\left\lbrack - k_{y}\left( b - b_{y}^{50} \right) \right\rbrack} \\
\end{matrix}
}{k_{y} = k \cdot \exp\left( \epsilon_{y,1}^{Sel} \right)
}{b_{y}^{50} = b^{50} \cdot exp(\epsilon_{y,2}^{Sel})}
```

where the parameters of the logistic form are allowed to vary over time.

In the context of semi-parametric selectivity, the following equation is
used:

``` math
{\begin{matrix}
{Sel}_{y,b}^{'} = \frac{1}{1 + \exp\left( - k\left\lbrack b - b^{50} \right\rbrack \right)}\exp\left( {\epsilon_{y,b}}^{Sel} \right) \\
\end{matrix}
}{{Sel}_{y,b} = \frac{{Sel}_{y,b}^{'}}{mean(\mathbf{Se}\mathbf{l}^{\mathbf{'}})}}
```

where deviations are placed about the parametric form and selectivity
values are mean standardized to aid with interpretability. Mean
standardization is applied only when semi-parametric deviations are
specified (process error models 3-5), or when non-parametric selectivity
is specified. For age-based selectivity, the mean is computed from a
single population and season reference ($`p = 1, \tau = 1`$) since the
underlying selectivity is invariant across these dimensions, and the
standardization is then applied identically across all populations and
seasons:

``` math
{Sel}_{r,y,b,s,j} = \exp\left(\log\left({Sel}_{r,y,b,s,j}\right) - \overline{\log\left(\mathbf{Sel}_{r,s,j}\right)}\right)
```

where $`\overline{\log\left(\mathbf{Sel}_{r,s,j}\right)}`$ is the mean
of log-selectivity across all years and bins for a given region, sex,
and fleet. For length-based selectivity, mean standardization is applied
directly to the selectivity-at-length values before conversion to the
age domain via the size-age transition matrix. Because the log-scale
non-parametric form already standardizes within each year over its own
window $`\mathcal{B}`$, it is excluded from this joint standardization
rather than being centered a second time. Further details on how
selectivity deviations arise can be found in the “Selectivity Process
Error” section of this document.

Finally, individual bins can be overridden with their own freely
estimated annual deviations (`_sel_bin_dev_bins`), independent of the
fleet’s functional form. For each named bin $`b^{}`$:

``` math
\text{Sel}_{r,y,b^{},s,j} = \exp\left( \epsilon_{r,y,b^{},s,j}^{\text{BinDev}} \right)
```

replacing whatever the functional form (and any standardization or
semi-parametric deviation) produced for that bin, while all remaining
bins keep their parametric shape. The override is applied last, after
every other transformation. The canonical use is a gear whose curve is
well described by a parametric form over most of its range but whose
youngest bin is governed by availability rather than by the gear (e.g.,
age-1 availability to a bottom trawl varying with year-class strength):
that bin becomes freely time-varying without abandoning the parametric
form elsewhere. The override deviations can carry their own iid or
random walk process error (`cont_tv_sel_bin_devs`), described in the
Selectivity Process Error section.

## Likelihoods

Currently, `SPoRC` incorporates data likelihood components for the
following data sources:

1.  region-aggregated fishery catches (summed across populations),
2.  population-specific fishery catches,
3.  region-aggregated fishery discards (summed across populations),
4.  population-specific fishery discards,
5.  region-aggregated fishery indices (summed across populations),
6.  population-specific fishery indices,
7.  region-aggregated fishery age compositions (summed across
    populations),
8.  population-specific fishery age compositions,
9.  region-aggregated fishery length compositions (summed across
    populations),
10. population-specific fishery length compositions,
11. region-aggregated discard age compositions (summed across
    populations),
12. population-specific discard age compositions,
13. region-aggregated discard length compositions (summed across
    populations),
14. population-specific discard length compositions,
15. region-aggregated survey indices (summed across populations),
16. population-specific survey indices,
17. region-aggregated survey age compositions (summed across
    populations),
18. population-specific survey age compositions,
19. region-aggregated survey length compositions (summed across
    populations),
20. population-specific survey length compositions, and
21. conventional tagging data.

Region-aggregated likelihoods compare observed data to predicted
quantities summed across all populations ($`\sum_p`$), while
population-specific likelihoods compare observed data to predicted
quantities for a single population $`p`$ directly. The total likelihood
(objective function) is the sum of the individual likelihood
contributions from these data sources along with priors and penalties,
where the objective function is minimized using a non-linear
optimization algorithm to estimate model parameters.

### Observation Likelihoods

#### Age-Disaggregated Observations

Catch, discards and the fishery and survey indices can each be fit at
age rather than aggregated with a composition alongside. Every age is
its own observation with its own standard deviation, and for the indices
its own catchability.

Each stream is stored over regions and sexes whatever a fleet reports,
and the fleet’s Type $`\mathcal{R}_{f}, \mathcal{S}_{f}`$ names which of
those margins it reports separately: a split margin is a single index, a
summed margin is the whole extent and the observation sits in its first
slot. Writing $`\mathcal{P} = \{1,\dots,n_{p}\}`$ for the populations,
the prediction for retained catch is

``` math
\text{CatchAA}_{r,y,\tau,a,\varsigma,f} =
\sum_{p \in \mathcal{P}} \sum_{r' \in \mathcal{R}_{f}(r)} \sum_{s \in \mathcal{S}_{f}(\varsigma)}
C_{p,r',y,\tau,a,s,f} \, \left[ w_{p,r',y,\tau,a,s,f} \right]^{u_{f}}
```

where $`\mathcal{R}_{f}(r) = \{r\}`$ when the fleet splits regions and
$`\{1,\dots,n_{r}\}`$ when it sums over them,
$`\mathcal{S}_{f}(\varsigma)`$ likewise for sexes, and
$`u_{f} \in \{0,1\}`$ selects abundance or biomass through the fleet’s
units. The population-specific form replaces $`\mathcal{P}`$ with the
single population being compared. The other two streams differ only in
the quantity summed:

``` math
\text{DiscardAA}_{r,y,\tau,a,\varsigma,f} = \sum \frac{D_{p,r',y,\tau,a,s,f}}{\text{dmr}_{r',y,\tau,f}}
\, \left[ w_{p,r',y,\tau,a,s,f} \right]^{u_{f}}
```

``` math
\text{SrvIdxAA}_{r,y,\tau,a,\varsigma,f} = \sum I_{p,r',y,\tau,a,s,f}
```

with the sums running over the same index sets. $`C`$ and $`D`$ are
retained catch and dead discards at age, $`\text{dmr}`$ the discard
mortality rate that raises the dead discards to the total discarded,
$`w`$ weight at age, $`N`$ numbers at age, $`S^{F}`$ and $`S^{R}`$
fishery selectivity and retention, and $`I`$ the survey-available
numbers at age. The indices carry no separate $`q`$: an index fit age by
age holds its age-specific catchability in selectivity, through the
`"nonparfree"` form.

A fleet’s observations are lognormal or normal. Writing $`\mu`$ for the
prediction and $`O`$ for the observation, the residual is
$`\varepsilon_{a} = \log O_{a} - \log \mu_{a}`$ under the lognormal and
$`\varepsilon_{a} = O_{a} - \mu_{a}`$ under the normal. The standard
deviation is whichever of the estimated parameter
$`\sigma_{a,\varsigma,f}`$ and the reported standard error
$`\varsigma^{obs}_{a}`$ the fleet’s error source names:

``` math
s_{a} = \sigma_{a,\varsigma,f}, \qquad \varsigma^{obs}_{a}, \qquad
\varsigma^{obs}_{a} + \sigma_{a,\varsigma,f}, \qquad
\sqrt{(\varsigma^{obs}_{a})^{2} + \sigma_{a,\varsigma,f}^{2}}
```

for `"none"`, `"data"`, `"est_additive"` and `"est_quadrature"`
respectively.

The standard deviations are coupled through integer key arrays over age,
sex and fleet, in which equal entries share a parameter and `NA`
excludes one. A single structure therefore gives one standard deviation
per age, one per age group, or one per fleet, and repeating entries
across sexes couples the sexes.

This is not the same statement as an aggregated observation with a
composition beside it. That factorization is exact for a multinomial
over Poisson counts, where the total is a genuine count sum, but the sum
of lognormals is not lognormal, so the two forms differ for the
lognormal used here. A fleet fits one or the other.

##### Correlation across ages

Ages within a cell may be independent, or correlated through one of
three structures, chosen per fleet. Under `"iid"` each age contributes
$`-\log \phi(\varepsilon_{a} ; 0, s_{a})`$ on its own. The other three
place the cell’s residual vector in one multivariate normal density,
$`\varepsilon \sim \text{MVN}(0, \Sigma)`$, and assign the whole density
to the first age present.

Under `"1dar1"` the correlation is a function of age distance, not of
position in the observed vector, so a fleet that skips ages is spaced by
the ages themselves:

``` math
\Sigma_{ij} = s_{a_i} \, s_{a_j} \, \rho_{r,\varsigma,f}^{\,|a_i - a_j|}
```

Where the observed ages are consecutive this is evaluated by the
autoregressive recursion, which is the same density at lower cost.

Under `"us"` the correlation is unstructured, built from
$`n_{a}(n_{a}-1)/2`$ unconstrained parameters that fill the strict lower
triangle of a matrix $`L`$ whose diagonal is one. Normalizing each row
of $`L`$ to unit length makes it a Cholesky factor, so

``` math
R = L L^{\top}, \qquad \Sigma_{ij} = s_{a_i} \, s_{a_j} \, R_{ij}
```

is a correlation matrix for any parameter values, and the principal
submatrix on whichever ages a cell observes is one too.

Under `"2dar1"` the correlation runs over ages and years jointly, with
the covariance the Kronecker product of an AR(1) over each:

``` math
\Sigma = \left( R^{\text{age}}(\rho_{r,\varsigma,f}) \otimes
R^{\text{year}}(\rho^{y}_{r,\varsigma,f}) \right) \odot \left( s \, s^{\top} \right),
\qquad R^{\text{age}}_{ij} = \rho^{\,|i-j|}
```

In all three the standard deviations enter as marginal standard
deviations, so $`s_{a}`$ means the same thing whichever structure a
fleet chooses and one key matrix serves all of them. This is worth
naming because the selectivity process error elsewhere in `SPoRC`
parameterizes its separable AR(1) by the conditional variance instead,
passing `dseparable` a scale of
$`\sigma / \sqrt{(1-\rho_{y}^{2})(1-\rho_{a}^{2})}`$; the two are
different parameterizations of the same family, not different families.

The correlations are indexed by region, sex and fleet, with a leading
population index for the population-specific streams, and are shared
through `rho_*_spec`.

This is defined on a complete grid, so a fleet’s observed ages and years
must form one. A cell with a single observed age has no correlation to
describe and falls back to independent.

The weight $`\lambda_{r,y,\tau,f}`$ is the one the aggregated stream
carries, applied after summing over ages and sexes within a cell.

#### Fishery Catches

Fishery catches can be fit using a lognormal likelihood. The
log-likelihood for region-aggregated observed catch,
$`\ell\left( \log\left( \text{ObsCatch}_{r,y,\tau,f} \right) \right)`$,
is defined as:

``` math
\begin{matrix}
\ell\left( \log\left( \text{ObsCatch}_{r,y,\tau,f} \right) \right) = \\
\end{matrix}
```

``` math
\lambda_{\text{ObsCatch}_{r,y,\tau,f}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsCatch}_{r,y,\tau,f}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsCatch}_{r,y,\tau,f} \right) - log\left( \text{Catch}_{r,y,\tau,f} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsCatch}_{r,y,\tau,f}}^{2}} \right)
```

Here, $`\lambda_{\text{ObsCatch}_{r,y,\tau,f}}`$ is the likelihood
weight, $`\text{ObsCatch}_{r,y,\tau,f}`$ is the observed catch,
$`\text{Catch}_{r,y,\tau,f}`$ is the predicted catch summed over
populations, and $`\sigma_{\text{ObsCatch}_{r,y,\tau,f}}^{2}`$ is the
variance of catch on the log scale.

Population-specific catch observations can additionally be fit using the
same lognormal form, comparing observed catch for a single population to
the predicted catch for that population without summing across
populations:

``` math
\begin{matrix}
\ell\left( \log\left( \text{ObsCatch}_{p,r,y,\tau,f} \right) \right) = \\
\end{matrix}
```

``` math
\lambda_{\text{ObsCatch}_{p,r,y,\tau,f}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsCatch}_{p,r,y,\tau,f}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsCatch}_{p,r,y,\tau,f} \right) - log\left( \text{Catch}_{p,r,y,\tau,f} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsCatch}_{p,r,y,\tau,f}}^{2}} \right)
```

where $`\text{Catch}_{p,r,y,\tau,f}`$ is the predicted catch for
population $`p`$ only.

#### Fishery Discards

Fishery discards are fit using the same lognormal likelihood form as
catches. The log-likelihood for region-aggregated observed discards is:

``` math
\begin{matrix}
\ell\left( \log\left( \text{ObsDiscard}_{r,y,\tau,f} \right) \right) = \\
\end{matrix}
```

``` math
\lambda_{\text{ObsDiscard}_{r,y,\tau,f}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsDiscard}_{r,y,\tau,f}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsDiscard}_{r,y,\tau,f} \right) - log\left( \text{Discard}_{r,y,\tau,f} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsDiscard}_{r,y,\tau,f}}^{2}} \right)
```

where $`\lambda_{\text{ObsDiscard}_{r,y,\tau,f}}`$ is the likelihood
weight, $`\text{ObsDiscard}_{r,y,\tau,f}`$ is the observed discard,
$`\text{Discard}_{r,y,\tau,f}`$ is the predicted discard summed over
populations, and $`\sigma_{\text{ObsDiscard}_{r,y,\tau,f}}^{2}`$ is the
variance of discards on the log scale. Population-specific discard
observations follow the same lognormal form with
$`\text{Discard}_{p,r,y,\tau,f}`$ for population $`p`$ only.

#### Fishery and Survey Indices

Fishery indices can also be fit assuming a lognormal likelihood. The
log-likelihood for region-aggregated observed fishery indices is:

``` math
\begin{matrix}
\ell\left( \log\left( \text{ObsFshIdx}_{r,y,\tau,f} \right) \right) = \\
\end{matrix}
```

``` math
\lambda_{\text{ObsFshIdx}_{r,y,\tau,f}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsFshIdx}_{r,y,\tau,f}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsFshIdx}_{r,y,\tau,f} \right) - log\left( \text{FshIdx}_{r,y,\tau,f} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsFshIdx}_{r,y,\tau,f}}^{2}} \right)
```

where $`\lambda_{\text{ObsFshIdx}_{r,y,\tau,f}}`$ controls the weight of
fishery indices to the objective function,
$`\text{ObsFshIdx}_{r,y,\tau,f}`$ represents the observed fishery
indices, $`\text{FshIdx}_{r,y,\tau,f}`$ is the predicted fishery index
summed across populations, and
$`\sigma_{\text{ObsFshIdx}_{r,y,\tau,f}}^{2}`$ denotes the variance of
the fishery index.

Population-specific fishery indices can additionally be fit, comparing
observed population-specific indices to the predicted index for that
population directly:

``` math
\begin{matrix}
\ell\left( \log\left( \text{ObsFshIdx}_{p,r,y,\tau,f} \right) \right) = \\
\end{matrix}
```

``` math
\lambda_{\text{ObsFshIdx}_{p,r,y,\tau,f}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsFshIdx}_{p,r,y,\tau,f}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsFshIdx}_{p,r,y,\tau,f} \right) - log\left( \text{FshIdx}_{p,r,y,\tau,f} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsFshIdx}_{p,r,y,\tau,f}}^{2}} \right)
```

where $`\text{FshIdx}_{p,r,y,\tau,f}`$ is the predicted fishery index
for population $`p`$ only.

Likewise, survey indices can be fit assuming a lognormal likelihood. The
log-likelihood for region-aggregated survey indices is:

``` math
\begin{matrix}
\ell\left( \log\left( \text{ObsSrvIdx}_{r,y,\tau,sf} \right) \right) = \\
\end{matrix}
```

``` math
\lambda_{\text{ObsSrvIdx}_{r,y,\tau,sf}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsSrvIdx}_{r,y,\tau,sf}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsSrvIdx}_{r,y,\tau,sf} \right) - log\left( \text{SrvIdx}_{r,y,\tau,sf} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsSrvIdx}_{r,y,\tau,sf}}^{2}} \right)
```

$`\lambda_{\text{ObsSrvIdx}_{r,y,\tau,sf}}`$ is the likelihood weight
applied to survey indices, $`\text{ObsSrvIdx}_{r,y,\tau,sf}`$ are the
observed survey indices, $`\text{SrvIdx}_{r,y,\tau,sf}`$ is the
predicted survey index summed across populations, and
$`\sigma_{\text{ObsSrvIdx}_{r,y,\tau,sf}}^{2}`$ indicates the variance
of the survey index.

Population-specific survey indices can additionally be fit, comparing
observed population-specific indices to the predicted index for that
population directly:

``` math
\begin{matrix}
\ell\left( \log\left( \text{ObsSrvIdx}_{p,r,y,\tau,sf} \right) \right) = \\
\end{matrix}
```

``` math
\lambda_{\text{ObsSrvIdx}_{p,r,y,\tau,sf}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsSrvIdx}_{p,r,y,\tau,sf}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsSrvIdx}_{p,r,y,\tau,sf} \right) - log\left( \text{SrvIdx}_{p,r,y,\tau,sf} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsSrvIdx}_{p,r,y,\tau,sf}}^{2}} \right)
```

where $`\text{SrvIdx}_{p,r,y,\tau,sf}`$ is the predicted survey index
for population $`p`$ only.

The lognormal above is the default, but each region-aggregated fishery
and survey index fleet can instead be assigned one of two alternative
error structures via `FishIdx_LikeType` / `SrvIdx_LikeType`. A normal
likelihood on the arithmetic scale treats the supplied standard errors
as arithmetic rather than log-scale:

``` math
\ell\left( \text{ObsIdx}_{r,y,\tau,j} \right) = \frac{1}{\sqrt{2\pi\sigma_{r,y,\tau,j}^{2}}\,}\exp\left( - \frac{\left\lbrack \text{ObsIdx}_{r,y,\tau,j} - \text{Idx}_{r,y,\tau,j} \right\rbrack^{2}}{2\sigma_{r,y,\tau,j}^{2}} \right)
```

A multivariate normal likelihood places a fleet’s whole observed series
in a single density with a fixed, user-supplied covariance
$`\mathbf{\Sigma}_{j}`$ across observations:

``` math
\ell\left( \mathbf{ObsIdx}_{j} \right) = \frac{1}{(2\pi)^{n/2}\left| \mathbf{\Sigma}_{j} \right|^{1/2}}\exp\left( - \frac{1}{2}\left( \mathbf{ObsIdx}_{j} - \mathbf{Idx}_{j} \right)^{T}\mathbf{\Sigma}_{j}^{- 1}\left( \mathbf{ObsIdx}_{j} - \mathbf{Idx}_{j} \right) \right)
```

where $`n`$ is the number of fitted observations for the fleet, ordered
as they appear when scanning the fleet’s use flags in array order. This
is the appropriate form when the survey itself provides a covariance
across years (e.g., a model-based index with estimated inter-annual
correlation), since a diagonal likelihood would treat correlated
residuals as independent information. The covariance is validated at
setup for symmetry and positive definiteness and is factorized once when
the AD tape is built. Note that the full multivariate density includes
the $`-\tfrac{n}{2}\log(2\pi)`$ constant, unlike TMB’s `MVNORM`
convention, so absolute likelihood values differ from implementations
that omit it even when the fits are identical. One-step-ahead (OSA)
residuals are available only for lognormal index fleets; normal and MVN
fleets are excluded from the OSA machinery.

#### Fishery and Survey Compositions

Several options for fitting composition data are available in `SPoRC`.
These include the multinomial, the Dirichlet-multinomial, and the
logistic-normal likelihoods. In the case of the multinomial likelihood,
the following expression is used:

``` math
{\begin{matrix}
\ell\left( \text{ObsCompositionData}_{r,y,\tau,j} \right) = \\
\end{matrix}
}{\lambda_{\text{ObsCompositionData}_{r,y,\tau,j}} \cdot \text{ISS}_{r,y,\tau,j} \cdot \prod_{b = 1}^{n_{b}}E_{r,y,\tau,b,j}^{O_{r,y,\tau,b,j}}}
```

where subscript $`j`$ is used to indicate a fishery or survey fleet and
the $`b`$ subscript generically indicates a bin number.
$`\lambda_{\text{ObsCompositionData}_{r,y,\tau,j}}`$ are likelihood
weights applied to composition data, $`ISS_{r,y,\tau,j}`$ is the input
sample size, $`E_{r,y,\tau,b,j}`$ denotes the expected composition
proportions, and $`O_{r,y,\tau,b,j}`$ are the observed composition
proportions.

In its negative log form the multinomial is implemented with a small
constant $`c`$ (`addtocomp`) guarding $`\log(0)`$ and an offset so a
perfect fit contributes zero:

``` math
-\ell = \text{ESS}\sum_{b}\left( O_{b} + c \cdot \mathbb{1}^{\text{const}} \right)\left\lbrack \log\left( O_{b} + c \right) - \log\left( E_{b} + c \right) \right\rbrack
```

The switch `comp_const_obs` ($`\mathbb{1}^{\text{const}}`$) controls
whether $`c`$ is also added to the observed proportions used as weights.
The default (`1`) is the long-standing `SPoRC` behavior; `0` adds the
constant only inside the logarithms, a convention several existing
assessments use, and matters when bridging to them, since the added
constant slightly reweights every bin.

If a Dirichlet-multinomial likelihood is assumed, the following
parameterization (linear) is used:

``` math
\begin{matrix}
\ell\left( \text{ObsCompositionData}_{r,y,\tau,j} \right) = \\
\lambda_{\text{ObsCompositionData}_{r,y,\tau,j}} \cdot \frac{\Gamma\left( \text{ISS}_{r,y,\tau,j} + 1 \right)}{\text{ISS}_{r,y,\tau,j} \cdot O_{r,y,\tau,b,j} + 1} \cdot \frac{\Gamma\left( \text{ISS}_{r,y,\tau,j}\theta_{r,j} \right)}{\text{Γ(}\text{ISS}_{r,y,\tau,j}\text{+}\text{ISS}_{r,y,\tau,j}\theta_{r,j}\text{)}} \cdot \\
\prod_{b = 1}^{n_{b}}\frac{\Gamma(\text{ISS}_{r,y,\tau,j}O_{r,y,\tau,b,j} + \text{ISS}_{r,y,\tau,j}\theta_{r,j}E_{r,y,\tau,b,j})}{\text{ISS}_{r,y,\tau,j}\theta_{r,j}E_{r,y,\tau,b,j}} \\
\end{matrix}
```

Here, $`\theta_{r,j}`$ is the overdispersion parameter of the
Dirichlet-multinomial that adjusts the input sample size. The effective
sample size ($`\text{ESS}_{r,y,\tau,j})`$ can then be derived as:

``` math
\begin{matrix}
\text{ESS}_{r,y,\tau,j} = \frac{1}{1 + \theta_{r,j}\ } + \text{ISS}_{r,y,\tau,j}\frac{\theta_{r,j}}{1 + \theta_{r,j}\ } \\
\end{matrix}
```

A multivariate logistic-normal likelihood can also be assumed, which is
given by:

``` math
{\begin{matrix}
\ell\left( \text{ObsCompositionData}_{r,y,\tau,j} \right) = \\
\end{matrix}
}{\frac{1}{(2\pi)^{\frac{B - 1}{2}}\left| \mathbf{\Sigma} \right|^{\frac{1}{2}}}\exp\left( - \frac{1}{2}\left\{ \mathbf{O}_{r,y,\tau,j}^{\mathbf{'}} - \mathbf{E}_{r,y,\tau,j}^{\mathbf{'}} \right\}^{T}\mathbf{\Sigma}^{- 1}\left\{ \mathbf{O}_{r,y,\tau,j}^{\mathbf{'}} - \mathbf{E}_{r,y,\tau,j}^{\mathbf{'}} \right\} \right)}
```

Both $`\mathbf{O}_{r,y,\tau,j}^{\mathbf{'}}`$ and
$`\mathbf{E}_{r,y,\tau,j}^{\mathbf{'}}`$ are $`(B - 1)`$ dimensional
vectors, while $`\mathbf{\Sigma}`$ is a $`(B - 1) \times (B - 1)`$
covariance matrix (see below for further details).
$`\mathbf{O}_{r,y,\tau,j}^{\mathbf{'}}`$ and
$`\mathbf{E}_{r,y,\tau,j}^{\mathbf{'}}`$ are derived via an additive
logistic function:

``` math
\begin{matrix}
O_{r,y,\tau,b,j}^{'} = \log\left( O_{r,y,\tau, - B,j} \right) - log(O_{r,y,\tau,B,j}) \\
E_{r,y,\tau,b,j}^{'} = \log\left( E_{r,y,\tau, - B,j} \right) - log(E_{r,y,\tau,B,j}) \\
\end{matrix}
```

where $`O_{r,y,\tau,b,j}^{'}`$ and $`E_{r,y,\tau,b,j}^{'}`$ are
transformed proportions using the last bin $`B`$ as the reference
category. Because the logarithm of zero is undefined, all untransformed
proportions must be strictly positive. If any observed proportion is
zero, both the observed and corresponding expected values are removed,
and the remaining proportions are renormalized to ensure that they sum
to one before applying the transformation. The covariance matrix of the
logistic-normal likelihood can be specified in various ways. In the
simplest case, the covariance matrix can be assumed to be independent
and identically distributed (iid):

``` math
\begin{matrix}
\mathbf{\Sigma =}\left( \mathbf{I}_{B} \cdot \theta^{2} \right)_{\mathbf{-}B} \\
\end{matrix}
```

where $`\mathbf{I}_{B}`$ is a $`B \times B`$ identity matrix and
$`\theta^{2}`$ is an estimated overdispersion parameter representing the
variance. The simple iid case can be further extended to incorporate a
one-dimensional lag-1 autoregressive structure:

``` math
\begin{matrix}
\mathbf{\Sigma =}\left( \mathbf{R}_{B} \cdot \frac{\theta^{2}}{1-\rho^2_B} \right)_{- B} \\
\left( \mathbf{R}_{B} \right)_{i,j}\mathbf{=}\rho_{B}^{|i - j|},\ \ i,j = \ 1,\ \cdots,\ B \\
\end{matrix}
```

Here, $`\mathbf{R}_{B}`$ is a $`B \times B`$ correlation matrix with a
lag-1 autoregressive structure, where $`\rho_{B}^{|i - j|}`$ defines the
correlation across bins. Lastly, if the model is specified to be
sex-structured and sex-composition data are utilized, a two-dimensional
autoregressive structure can be specified:

``` math
\begin{matrix}
\mathbf{\Sigma =}\left( {\mathbf{R}_{S}\mathbf{\ \bigotimes\ R}}_{C} \cdot \frac{\theta^{2}}{(1-\rho^2_S)(1-\rho^2_C)} \right)_{- B} \\
\end{matrix}
```

``` math
\left( \mathbf{R}_{C} \right)_{i,j}\mathbf{=}\rho_{C}^{|i - j|},\ \ i,j = \ 1,\ \cdots,\ C
```

``` math
\left( \mathbf{R}_{S} \right)_{i,j}\mathbf{=}\left\{ \begin{matrix}
1,\ \ if i = j \\
\rho_{s},\ \ if i \neq j \\
\end{matrix} \right.\ ,\ \ i,j = 1,\ldots,n_{s}
```

$`\mathbf{R}_{S}`$ is a constant correlation matrix dimensioned by
$`n_{s} \times n_{s}`$ for sexes, with off-diagonal elements
$`\rho_{s}`$ controlling the correlation of age/length categories across
sexes, while $`\mathbf{R}_{C}`$ is a $`n_{c} \times n_{c}`$ lag-1
autoregressive correlation structure, where $`\rho_{C}^{|i - j|}`$
defines the correlation across age/length categories.
$`\mathbf{\bigotimes}`$ denotes the Kronecker product.

All three composition likelihood forms (multinomial,
Dirichlet-multinomial, logistic-normal) can be applied to retained
fishery, discarded fishery, and survey composition data, as well as to
both region-aggregated and population-specific variants. For
region-aggregated compositions, expected values $`E_{r,y,\tau,b,j}`$ are
derived from catch-at-age or survey index-at-age quantities summed
across populations ($`\sum_p`$). For population-specific compositions,
expected values $`E_{p,r,y,\tau,b,j}`$ are derived from the quantities
for a single population $`p`$ directly. For discard compositions,
expected values are derived from discarded catch-at-age
($`D_{p,r,y,\tau,a,s,f}^{a}`$) or discarded catch-at-length quantities
analogously. Each likelihood form and covariance structure described
above applies identically across all composition data types;
population-specific likelihoods additionally carry separate
overdispersion ($`\theta_{p,r,j}`$) and correlation parameters
($`\rho_{p,r,j}`$) estimated independently from their region-aggregated
counterparts.

##### Structuring Compositions and Ageing Error

Related to the use of composition data likelihoods, composition data can
be structured differently depending on model assumptions and data
constraints. In particular, three options are available to fit to
composition data:

1.  ‘Aggregated’ compositions across regions and sexes,

2.  ‘Split’ compositions for each region and sex (i.e., no implicit
    information about sex-ratios), and

3.  ‘Joint’ compositions across sexes (i.e., implicit information is
    provided about sex-ratios).

The expected compositions (i.e., catch-at-age, catch-at-length, survey
catch-at-age, survey catch-at-length) when specified as ‘aggregated’ are
derived with the following:

``` math
{\begin{matrix}
E_{y,\tau,b}^{'} = \frac{\sum_{r = 1}^{n_{r}}{\sum_{s = 1}^{n_{s}}E_{r,y,\tau,b,s}^{'}}}{n_{s} \cdot n_{r}} \\
\end{matrix}
}{\mathbf{E}_{y,\tau} = \mathbf{E}_{y,\tau}^{\mathbf{'}}\mathbf{\Theta}_{y}}
```

where compositions are summed across regions and sexes and normalized to
sum to one. Ageing error ($`\mathbf{\Theta}_{y}`$) can then be applied
using standard matrix multiplication. Expected compositions that are
specified as ‘Split’ by sexes and regions are computed as:

``` math
\begin{matrix}
E_{y,\tau,b,s}^{'} = \frac{E_{r,y,\tau,b,s}^{'}}{\sum_{b = 1}^{n_{B}}E_{r,y,\tau,b,s}^{'}} \\
\end{matrix}
```

``` math
\mathbf{E}_{y,\tau,s} = \mathbf{E}_{y,\tau,s}^{\mathbf{'}}\mathbf{\Theta}_{y}
```

Here, expected compositions sum to one within a given region and sex
combination and ageing error is similarly applied via matrix
multiplication. In the case where expected compositions are specified as
‘Joint’, they are calculated as:

``` math
\begin{matrix}
E_{y,\tau,b,s}^{'} = \frac{E_{r,y,\tau,b,s}^{'}}{\sum_{s = 1}^{n_{s}}{\sum_{b = 1}^{n_{B}}E_{r,y,\tau,b,s}^{'}}} \\
\end{matrix}
```

``` math
\mathbf{E}_{y,\tau} = \left( \mathbf{E}_{y,\tau} \right)^{T}(\mathbf{I}_{s}\mathbf{\ \bigotimes\ \Theta}_{y})
```

where the expected compositions sum to one jointly across bins and
sexes, thus preserving implicit sex-ratio information. Ageing error is
then applied by taking the Kronecker product of a
$`n_{s}\ \times\ n_{s}`$ identity matrix with the ageing error matrix,
followed by matrix multiplication. These three structuring options apply
identically to retained fishery, discarded fishery, and survey
composition likelihoods, as well as to both region-aggregated and
population-specific variants.

Ageing error is fleet specific, $`\mathbf{\Theta}_{y,f}`$, supplied
through `AgeingError_fish` and `AgeingError_srv` to accomodate different
reading mehtods. Every fleet’s matrix maps onto the same set of observed
age bins, because the observed composition arrays carry one age
dimension shared across fleets. The length analog $`\mathbf{\Lambda}`$
(`LenBinMap`) is applied in exactly the same position for length
compositions, and is shared across fleets.

Additionally, every composition stream can be restricted to a subset of
observed bins (the `_bins` arguments). Both the observed and expected
compositions are subset to the named bins and renormalized within them:

``` math
E_{b}^{\text{fit}} = \frac{E_{b}}{\sum_{b' \in \mathcal{B}}E_{b'}},\qquad O_{b}^{\text{fit}} = \frac{O_{b}}{\sum_{b' \in \mathcal{B}}O_{b'}},\qquad b \in \mathcal{B}
```

where $`\mathcal{B}`$ is the fitted bin set, indexed on observed bins
(i.e., after $`\mathbf{\Theta}_{y,f}`$ or $`\mathbf{\Lambda}`$ has
mapped model bins onto observed ones). Bins outside $`\mathcal{B}`$ are
left out of the likelihood entirely rather than being forced to be
explained, which is appropriate for a gear that only resolves part of
the range (e.g., a fishery that never encounters the youngest ages,
whose zeros would otherwise carry information).

The restriction and the mapping are distinct operations applied in that
order. $`\mathbf{\Theta}_{y,f}`$ and $`\mathbf{\Lambda}`$ redistribute
mass from model bins onto observed bins, conserving it wherever a row
sums to one, which is every row of an ageing error matrix or a length
bin collapse. Restricting to $`\mathcal{B}`$ discards the mass outside
$`\mathcal{B}`$ and renormalizes what remains. A row that sums to zero
is the one case where a mapping also discards, which is what the shifted
identity $`\mathbf{I}[, 2{:}n_{A}]`$ does to drop a model age the
observations never resolve. Neither substitutes for the other: no subset
of the model’s bins reproduces a many-to-one collapse, and a mapping
that zeroed an observed bin would leave it with a structural-zero
expectation the composition likelihood cannot fit.

For ‘Joint’ compositions the restriction is taken on each sex’s block of
the $`\left\lbrack b \times s \right\rbrack`$ stack, so the implicit sex
ratio those compositions carry becomes the ratio within $`\mathcal{B}`$.
Where the logistic-normal families supply a correlation over bins, the
covariance is formed across all observed bins and then subset to
$`\mathcal{B}`$, so a gap in $`\mathcal{B}`$ still contributes to the
lag between the bins either side of it.

##### Conditional Age-at-Length

A conditional age-at-length observation is the age composition of the
fish aged from one length bin. Its expectation is the row of the joint
array for that bin, normalized across ages, then passed through the
ageing error matrix:

``` math
E_{p,r,y,\tau,a,s,f}^{\,\text{caal}(l)} = \dfrac{C_{p,r,y,\tau,l,a,s,f}^{la}}{\sum_{a'}C_{p,r,y,\tau,l,a',s,f}^{la}}
```

which is $`P(a \mid l)`$ under the model, since the joint array is
$`P(l \mid a)`$ times the numbers at age. Each length bin is then fit as
its own composition with the multinomial or Dirichlet-multinomial above,
with the input sample size the number of fish aged from that bin
(`ISS_Fish_caal`, `ISS_Srv_caal`) and the weight `Wt_Fish_caal`,
`Wt_Srv_caal` multiplying it:

``` math
-\ell^{\text{caal}} = \sum_{l}\lambda_{l}\,\text{ISS}_{l}\sum_{a}\left( O_{a}^{(l)} + c\,\mathbb{1}^{\text{const}} \right)\left\lbrack \log\left( O_{a}^{(l)} + c \right) - \log\left( E_{a}^{(l)} + c \right) \right\rbrack
```

with $`\theta`$ shared across the length bins of a fleet under the
Dirichlet-multinomial, since the bins come from one length-stratified
sample. The logistic-normal families are not available: a bin’s age
sample is small and mostly zeros, which the additive log-ratio transform
cannot handle. Fit together with the marginal length compositions, the
lengths carry the abundance signal and the conditional rows carry
$`P(a \mid l)`$, which is the information on growth; fitting marginal
ages as well would count the same fish twice.

#### Tagging

`SPoRC` currently allows for various tagging likelihoods, ranging from
the Poisson, Negative Binomial, multinomial, and Dirichlet-multinomial
likelihood. Additionally, `SPoRC` also allows for both release- and
recapture-conditioned dynamics (McGarvey and Feenstra, 2002). The
Poisson tag likelihood is given by:

``` math
\begin{matrix}
\ell\left( {ObsRecap}_{p,r,y,\tau,a,s,f}^{k} \right) = \lambda_{\text{Tagging}}\frac{\exp\left( - {Recap}_{p,r,y,\tau,a,s,f}^{k} \right)\left( {Recap}_{p,r,y,\tau,a,s,f}^{k} \right)^{{ObsRecap}_{p,r,y,\tau,a,s,f}^{k}}}{{ObsRecap}_{p,r,y,\tau,a,s,f}^{k}!} \\
\end{matrix}
```

where $`{ObsRecap}_{p,r,y,\tau,a,s,f}^{k}`$ are the observed tag
recaptures and $`\lambda_{\text{Tagging}}`$ is the likelihood weight
applied to tagging data. In the case where the Negative Binomial is
invoked, the following expression is used:

``` math
\begin{matrix}
\ell\left( {ObsRecap}_{p,r,y,\tau,a,s,f}^{k} \right) = \\
\lambda_{\text{Tagging}}\frac{\Gamma\left( {ObsRecap}_{p,r,y,\tau,a,s,f}^{k} + \eta \right)}{\Gamma(\eta)\Gamma\left( {ObsRecap}_{p,r,y,\tau,a,s,f}^{k} + 1 \right)}\left( \frac{\eta}{{Recap}_{p,r,y,\tau,a,s,f}^{k} + \eta} \right)^{\eta}\left( \frac{{Recap}_{p,r,y,\tau,a,s,f}^{k}}{{Recap}_{p,r,y,\tau,a,s,f}^{k} + \eta} \right)^{{ObsRecap}_{p,r,y,\tau,a,s,f}^{k}} \\
\end{matrix}
```

Here, $`\eta`$ represents the estimated overdispersion parameter for
tagging data.

Under release conditioned dynamics, both recaptured and non-recaptured
states are fit to. Proportions of observed
($`{PObsRecap}_{p,r,y,\tau,a,s,f}^{k}`$) and expected recaptured
($`{PRecap}_{p,r,y,\tau,a,s,f}^{k})`$ individuals are given by:

``` math
\begin{matrix}
{PObsRecap}_{p,r,y,\tau,a,s,f}^{k} = \frac{{ObsRecap}_{p,r,y,\tau,a,s,f}^{k}}{{InitTag}^{k}} \\
\end{matrix}
```

``` math
{PRecap}_{p,r,y,\tau,a,s,f}^{k} = \frac{{Recap}_{p,r,y,\tau,a,s,f}^{k}}{{InitTag}^{k}}
```

$`{InitTag}^{k}`$ denotes the total tags released for a given tag cohort
(combination of release region, year, and season). Non-recaptured states
can then be written as:

``` math
\begin{matrix}
{PObsNonRecap}^{k} = 1 - \sum_{p}^{}{\sum_{r}^{}{\sum_{y}^{}{\sum_{\tau}^{}{\sum_{a}^{}{\sum_{s}^{}{\sum_{f}^{}{PObsRecap}_{p,r,y,\tau,a,s,f}^{k}}}}}}} \\
\end{matrix}
```

``` math
{PNonRecap}^{k} = 1 - \sum_{p}^{}{\sum_{r}^{}{\sum_{y}^{}{\sum_{\tau}^{}{\sum_{a}^{}{\sum_{s}^{}{\sum_{f}^{}{PRecap}_{p,r,y,\tau,a,s,f}^{k}}}}}}}
```

where $`{PObsNonRecap}^{k}`$ and $`{PNonRecap}^{k}`$ are the observed
and expected non-recaptured states, respectively. These states are then
combined into a single vector of observed and expected values:

``` math
\begin{matrix}
\mathbf{O}_{Tagging}^{k} = \left\{ \mathbf{PObsReca}\mathbf{p}^{k},{PObsNonRecap}^{k} \right\} \\
\end{matrix}
```

``` math
\mathbf{E}_{Tagging}^{k} = \{\mathbf{PReca}\mathbf{p}^{k},{PNonRecap}^{k}\}
```

If a Multinomial likelihood is assumed for release conditioned dynamics,
this is given by:

``` math
\begin{matrix}
\ell\left( {ObsRecap}^{k} \right) = {\lambda_{\text{Tagging}}InitTag}^{k}\prod_{i}^{}\left( E_{i,Tagging}^{k} \right)^{O_{i,Tagging}^{k}} \\
\end{matrix}
```

Here, the subscript $`i`$ is used to generically denote a given element.
If a Dirichlet-multinomial with released-condition dynamics was assumed,
the tagging likelihood would be written as:

``` math
\begin{matrix}
\ell\left( {ObsRecap}^{k} \right) = \\
\lambda_{\text{Tagging}} \cdot \frac{\Gamma\left( {InitTag}^{k}\  + 1 \right)}{{InitTag}^{k}\  \cdot O_{i,Tagging}^{k} + 1} \cdot \frac{\Gamma\left( {InitTag}^{k}\ \eta \right)}{\text{Γ(}{InitTag}^{k}\text{+}{InitTag}^{k}\eta\text{)}} \cdot \\
\prod_{i}^{}\frac{\Gamma({InitTag}^{k} \cdot O_{i,Tagging}^{k} + {InitTag}^{k} \cdot \eta \cdot E_{i,Tagging}^{k})}{{InitTag}^{k} \cdot \eta \cdot E_{i,Tagging}^{k}} \\
\end{matrix}
```

The $`\eta`$ parameter in the Dirichlet-multinomial likelihood
represents the overdispersion parameter for tagging data.

Under recapture-conditioned dynamics, tag shedding, tag induced
mortality, and tag reporting rates are assumed to be spatially-invariant
and do not need to be estimated, given that these terms cancel out in
the denominator (McGarvey and Feenstra, 2002). Unlike
release-conditioned dynamics, assuming recaptured-conditioned processes
does not require fitting to non-recaptured states. Thus, the observed
and expected recaptured proportions can be written as:

``` math
{\begin{matrix}
{PObsRecap}_{p,r,y,\tau,a,s,f}^{k} = \frac{{ObsRecap}_{p,r,y,\tau,a,s,f}^{k}}{\sum_{p}^{}{\sum_{r}^{}{\sum_{a}^{}{\sum_{s}^{}{ObsRecap}_{p,r,y,\tau,a,s,f}^{k}}}}} \\
\end{matrix}
}{{PRecap}_{p,r,y,\tau,a,s,f}^{k} = \frac{{Recap}_{p,r,y,\tau,a,s,f}^{k}}{\sum_{p}^{}{\sum_{r}^{}{\sum_{a}^{}{\sum_{s}^{}{Recap}_{p,r,y,\tau,a,s,f}^{k}}}}}}
```

where recapture probabilities are normalized by the total number of
recaptures across populations, regions, ages, and sexes in a given year
and season.

### Parameter Priors and Process Error Penalties

#### Parameter Priors

Considering the complexity of integrated population models, several
priors can be specified to help inform the estimation of parameters by
providing additional knowledge. Priors can currently be specified for
natural mortality, fishery and survey catchability, fishery and survey
selectivity, steepness, recruitment population scale ($`R_0`$) and
proportions, stray rates, movement rates, and tag reporting rates.

##### Natural Mortality

In the case of natural mortality, a lognormal prior is utilized:

``` math
\begin{matrix}
P\left( \log\left( \text{NatMort}_{p,r,y,a,s} \right) \right) = \frac{1}{\sqrt{2\pi\sigma_{p,r,y,a,s}^{2\ (Natmort)}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{NatMort}_{p,r,y,a,s} \right) - log\left( \mu_{p,r,y,a,s}^{\text{(NatMort)}} \right) \right\rbrack^{2}}{2\sigma_{p,r,y,a,s}^{2\ (Natmort)}} \right) \\
\end{matrix}
```

where the variance of the prior is given by
$`\sigma_{p,r,y,a,s}^{2\ (Natmort)}`$, and
$`\mu_{p,r,y,a,s}^{\text{(NatMort)}}`$ denotes the prior mean.

##### Fishery and Survey Catchability

For fishery and survey catchability, a lognormal prior can also be
specified:

``` math
\begin{matrix}
P\left( \log\text{(}q_{r,y,j}\text{)} \right) = \frac{1}{\sqrt{2\pi\sigma_{r,y,j}^{2\ (qPrior)}}\,}\exp\left( - \frac{\left\lbrack \log\left( q_{r,y,j} \right) - log\left( \mu_{r,y,j}^{\text{(qPrior)}} \right) \right\rbrack^{2}}{2\sigma_{r,y,j}^{2\ (qPrior)}} \right) \\
\end{matrix}
```

where $`j`$ here represents either a fishery or survey fleet,
$`\sigma_{r,y,j}^{2\ (qPrior)}`$ is the variance of the prior, and
$`\mu_{r,y,j}^{\text{(qPrior)}}`$ indicates the prior mean for
catchability.

##### Fishery and Survey Selectivity

In general, selectivity priors can be utilized to serve as regularizing
priors to facilitate stable parameter estimation (Monnahan, 2024). These
priors are assumed to be lognormal and are applied to the selectivity
parameters themselves:

``` math
\begin{matrix}
P\left( \log\left( \theta_{r,p,y,s,j} \right) \right) = \frac{1}{\sqrt{2\pi{\sigma^{2}}_{r,p,y,s,j}^{\text{(Sel)}}}\,}\exp\left( - \frac{\left\lbrack \ln\left( \theta_{r,p,y,s,j} \right) - ln\left( \mu_{r,p,y,s,j}^{\text{(Sel)}} \right) \right\rbrack^{2}}{2{\sigma^{2}}_{r,p,y,s,j}^{\text{(Sel)}}} \right) \\
\end{matrix}
```

where $`\theta_{r,p,y,s,j}`$ is a selectivity parameter for a given
functional form specified, $`{\sigma^{2}}_{r,p,y,s,j}^{\text{(Sel)}}`$
is the prior variance, and $`\mu_{r,p,y,s,j}^{\text{(Sel)}}`$ is the
prior mean for the specific selectivity parameter.

In addition to lognormal priors on individual parameters, a centering
penalty can be applied jointly to a named set of selectivity
fixed-effect parameters (`Use__selex_penalty` with a `_selex_penalty`
table). For a set $`\mathcal{P}`$ of parameters held on the log scale,
the penalty is:

``` math
P_{\text{SelCenter}} = w\left\lbrack \log\left( \frac{1}{|\mathcal{P}|}\sum_{i \in \mathcal{P}}\exp\left( \theta_{i} \right) \right) \right\rbrack^{2}
```

i.e., the squared log of the set’s mean selectivity on the natural
scale, pushing the average selectivity of the set toward one. This
resolves the scale non-identifiability of non-parametric log-scale
selectivity: such a curve is only identified up to a scalar once
catchability or fishing mortality is free to absorb its level, and the
centering penalty pins that scalar softly rather than fixing a bin
outright. Because the expression averages on the natural scale, it is
meaningful only for parameter sets held on the log scale (the
`"nonparlog"` form); logit-scale sets would not average to anything
interpretable as selectivity.

##### Steepness

If a stock-recruit relationship is assumed (Beverton-Holt or Ricker),
priors for steepness can be specified. Currently, a scaled beta prior
(bounded between 0.2 and 1 by default) can be invoked:

``` math
{\begin{matrix}
a_{p,r}^{(h)} = \left( \frac{\mu_{p,r}^{\text{(}\text{h}\text{)}} - 0.2}{1 - 0.2} \right)\left( \frac{\sigma_{p,r}^{(h)}}{1 - 0.2} \right)^{2} \\
\end{matrix}
}{b_{p,r}^{(h)} = \left\lbrack 1 - \ \left( \frac{\mu_{p,r}^{\text{(}\text{h}\text{)}} - 0.2}{1 - 0.2} \right) \right\rbrack\left\lbrack \left( \frac{\sigma_{p,r}^{(h)}}{1 - 0.2} \right)^{2} \right\rbrack
}{P\left( h_{p,r} \right) = \frac{\Gamma\left( a_{p,r}^{(h)} \right)\Gamma\left( b_{p,r}^{(h)} \right)}{\Gamma\left( a_{p,r}^{(h)} + b_{p,r}^{(h)} \right)}h_{p,r}^{a - 1}\left( 1 - h_{p,r} \right)^{b - 1}}
```

Here, $`a_{p,r}^{(h)}`$ and $`b_{p,r}^{(h)}`$ are parameters of the beta
distribution, $`\mu_{p,r}^{\text{(}\text{h}\text{)}}`$ is the prior mean
steepness for a given population and region (bounded between 0.2 and 1)
while $`\sigma_{p,r}^{(h)}`$ is the standard deviation for these priors.

The support of the beta can be changed from the default $`(0.2, 1)`$ via
optional `lb` and `ub` columns of the prior table, in which case every
$`0.2`$ and $`1`$ above is replaced by the supplied bounds. This exists
because a beta placed on a different interval is a genuinely different
function of $`h`$ rather than the same one shifted: a beta on $`(0, 1)`$
carries a $`\log(h)`$ term where the rescaled default carries
$`\log(h - 0.2)`$, and no choice of shape parameters reconciles the two.
Matching a bridged assessment’s steepness prior therefore requires
matching its support, not just its mean and standard deviation.

##### Recruitment Proportions

Regional recruitment is derived by apportioning a global recruitment
parameter using regional recruitment proportions for each population
(i.e., $`\mu_p^{\text{Rec}} \cdot \zeta_{p,r}`$). Here, $`\zeta_{p,r}`$
is derived via a multinomial logit transformation and Dirichlet priors
can be used to help constrain estimation:

``` math
\begin{matrix}
P\left( \mathbf{\zeta}_p \right) = \frac{\Gamma\left( \sum_{r}^{n_{r}}\alpha_{p,r} \right)}{\prod_{r}^{n_{r}}{\Gamma(\alpha_{p,r})}}\prod_{r = 1}^{n_{r}}\zeta_{p,r}^{\alpha_{p,r} - 1} \\
\end{matrix}
```

$`\mathbf{\zeta}_p = \{\zeta_{p,1},\zeta_{p,2},\ldots,\zeta_{p,n_{r}}\}`$
are the estimated recruitment proportions across regions for population
$`p`$, and $`\alpha_{p,r}`$ is the concentration parameter governing the
spread of the Dirichlet distribution. Similarly, seasonal recruitment
proportions $`\chi_{p,\tau}`$ can be constrained with Dirichlet priors
when estimated. When $`RecLag = 0`$ and $`\tau^{spawn} > 1`$,
$`\chi_{p,\tau}`$ is instead parameterized via a multinomial logit
restricted to seasons $`\tau^{spawn},\ldots,n_\tau`$ (seasons before
$`\tau^{spawn}`$ are fixed at exactly zero rather than estimated, per
the timing constraint described under Age-0 Recruitment above), and any
Dirichlet prior is evaluated only over that same restricted support.

##### R0

A lognormal prior can be placed on $`R_0`$ for any population:

``` math
P\left(\ln R_{0,p}\right) = \frac{1}{\sigma_p^{(R_0)}\sqrt{2\pi}}\exp\left(-\frac{\left(\ln R_{0,p} - \ln \mu_p^{(R_0)}\right)^2}{2\left(\sigma_p^{(R_0)}\right)^2}\right)
```

where $`\mu_p^{(R_0)}`$ is the prior mean on the natural scale and
$`\sigma_p^{(R_0)}`$ is the standard deviation on the log scale.

##### Stray Rates

When stray rates are estimated ($`n_p > 1`$ and
`use_fixed_stray_rate = 0`), a standard beta prior can be applied to
regularize estimation. The prior is parameterized via method-of-moments
in terms of a mean and standard deviation:

``` math
\kappa_{p,b}^{(\phi)} = \frac{\mu_{p,b}^{(\phi)}\left(1 - \mu_{p,b}^{(\phi)}\right)}{\left(\sigma_{p,b}^{(\phi)}\right)^{2}} - 1
```

``` math
a_{p,b}^{(\phi)} = \mu_{p,b}^{(\phi)} \cdot \kappa_{p,b}^{(\phi)}
```

``` math
b_{p,b}^{(\phi)} = \left(1 - \mu_{p,b}^{(\phi)}\right) \cdot \kappa_{p,b}^{(\phi)}
```

where $`\kappa_{p,b}^{(\phi)}`$ is the concentration parameter. The
stray rate is numerically stabilized by squishing the logistic transform
away from the boundaries:

``` math
\tilde{\phi}_{p,b} = \epsilon + (1 - 2\epsilon) \cdot \text{logistic}\left(\phi_{p,b}^{}\right)
```

where $`\phi_{p,b}^{}`$ is the logit-scale parameter and $`\epsilon`$ is
a small constant (e.g. $`10^{-4}`$) ensuring
$`\tilde{\phi}_{p,b} \in (\epsilon, 1-\epsilon)`$. The prior is then:

``` math
P\left(\tilde{\phi}_{p,b}\right) = \frac{\Gamma\left(a_{p,b}^{(\phi)} + b_{p,b}^{(\phi)}\right)}{\Gamma\left(a_{p,b}^{(\phi)}\right)\Gamma\left(b_{p,b}^{(\phi)}\right)}\tilde{\phi}_{p,b}^{a_{p,b}^{(\phi)} - 1}\left(1 - \tilde{\phi}_{p,b}\right)^{b_{p,b}^{(\phi)} - 1}
```

where $`\sigma_{p,b}^{(\phi)}`$ is the literal standard deviation of the
Beta distribution and must satisfy
$`\sigma_{p,b}^{(\phi)} < \sqrt{\mu_{p,b}^{(\phi)}\left(1 - \mu_{p,b}^{(\phi)}\right)}`$
to ensure $`\kappa_{p,b}^{(\phi)} > 0`$. Because stray rates are
generally not identifiable from fisheries data alone, this prior serves
primarily as a regularizing constraint rather than an informative prior,
and tight values of $`\sigma_{p,b}^{(\phi)}`$ are recommended. Note that
when $`\sigma_{p,b}^{(\phi)}`$ is large relative to
$`\mu_{p,b}^{(\phi)}(1 - \mu_{p,b}^{(\phi)})`$, $`a`$ and $`b`$ approach
zero and the Beta density becomes U-shaped, placing mass near 0 and 1.
In this regime numerical instability can occur during optimization,
which is why $`\tilde{\phi}_{p,b}`$ is squished away from the boundaries
via the $`\epsilon`$ transformation.

##### Movement

Likewise, priors on movement values can be assumed to arise from a
Dirichlet process:

``` math
\begin{matrix}
P\left( \mathbf{M}_{p,.,k,y,\tau,a,s} \right) = \frac{\Gamma\left( \sum_{r}^{n_{r}}c_{p,r,k,y,\tau,a,s} \right)}{\prod_{r}^{n_{r}}{\Gamma\left( c_{p,r,k,y,\tau,a,s} \right)}}\prod_{r = 1}^{n_{r}}M_{p,r,k,y,\tau,a,s}^{c_{p,r,k,y,\tau,a,s} - 1} \\
\end{matrix}
```

where $`r`$ is the origin region, $`k`$ is the destination, and
$`c_{p,r,k,y,\tau,a,s}`$ are the concentration parameters that control
the Dirichlet distribution.

###### What the prior is evaluated on

For unstructured movement (`move_type = 0`) the $`\mathbf{M}`$ above are
the estimated transition fractions themselves, and the prior applies to
them directly. For CTMC movement (`move_type = 1`) the model estimates a
generator $`\dot{\mathbf{Q}}`$ rather than fractions, and the prior is
evaluated on the annual fractions obtained from it,

``` math
\mathbf{m}_{p,r,\cdot,y,\tau,a,s} = \left\lbrack \exp\left( \dot{\mathbf{Q}}_{p,y,\tau,a,s} \right) \right\rbrack_{r,\cdot}
```

i.e. at $`\Delta t = 1`$, rather than on the seasonal fractions
$`\exp(\dot{\mathbf{Q}}\Delta t)`$ stored in $`\mathbf{M}`$. The
distinction matters once `ctmc_scale_by_seasdur = 1`: a season’s
movement matrix approaches the identity as the season shortens, so a
fixed $`\mathbf{c}`$ would silently become a far stronger constraint as
$`n_{\tau}`$ grows. On a three-region example the same $`c = 3`$ prior
costs 1.04 nLL units at $`n_{\tau} = 1`$ but 9.91 at $`n_{\tau} = 12`$;
evaluated annually it is 1.04 in both cases. Fixing the exponent at one
year also matches how such priors are elicited, as a belief about the
fraction of fish moving per year, independent of how the year is
partitioned. Under the legacy `ctmc_scale_by_seasdur = 0` the two
coincide.

Because $`\exp(\dot{\mathbf{Q}})`$ is a deterministic function of the
estimated CTMC parameters, the prior still informs $`\dot{\mathbf{Q}}`$
under every `move_timing`. Note, however, that under `move_timing = 2`
the seasonal operator is
$`\exp\left( \dot{\mathbf{Q}}^{T}\Delta\tau - \text{diag}(\mathbf{Z}) \right)`$,
which does not factor into movement and survival: $`\mathbf{M}`$ is then
a reported diagnostic rather than a term in the dynamics, and the prior
acts as a reparameterized prior on the generator.

###### Choosing the concentration parameters

Writing $`\mathbf{c}`$ for the concentration vector of one origin region
and $`n_{r}`$ for the number of regions:

``` math
\text{mode}_{k} = \frac{c_{k} - 1}{\sum_{j} c_{j} - n_{r}}\ \ (c_{k} > 1), \qquad \text{mean}_{k} = \frac{c_{k}}{\sum_{j} c_{j}}, \qquad \kappa = \sum_{j} c_{j} - n_{r}
```

where $`\kappa`$ behaves as the number of prior pseudo-observations, so
it is the natural handle on prior strength. Two cases are worth
separating:

- $`\mathbf{c} = \mathbf{1}`$ gives a uniform density over the simplex:
  every movement vector is equally likely, the contribution to the
  objective is constant, and the gradient is exactly zero. This is the
  uninformative choice.
- $`\mathbf{c} = c\mathbf{1}`$ with $`c > 1`$ is symmetric with its mode
  at $`1/n_{r}`$ in every region, a prior centered on equal movement
  among all regions, which is informative and, for a residency-dominated
  stock, a strong assumption rather than a vague one.

To target a particular movement vector with a chosen strength, set
$`c_{k} = 1 + \kappa\,\text{mode}_{k}`$.

##### Tag Reporting Rates

Two types of priors can be specified for tag reporting rates. In
particular, a symmetric beta distribution is applied:

``` math
\begin{matrix}
P\left( \beta_{r,y,f} \right) = \left( \beta_{r,y,f} + 1e - 4 \right)^{\sigma_{r,y,f}^{(\beta)}}\left( 1 - \beta_{r,y,f} + 1e - 4 \right)^{\sigma_{r,y,f}^{(\beta)}} \\
\end{matrix}
```

Here, $`\sigma_{r,y,f}^{(\beta)}`$ determines the scale of the tag
reporting parameter and determines how strongly to penalize estimates
when they approach the bounds of $`\lbrack 0,1\rbrack`$. Smaller values
of $`\sigma_{r,y,f}^{(\beta)}`$ result in larger penalties, and vice
versa.

Tag reporting rate priors can also be specified as a standard beta
distribution, parameterized via method-of-moments in terms of a mean and
standard deviation:

``` math
\kappa_{r,f}^{(\beta)} = \frac{\mu_{r,y,f}^{(\beta)}\left(1 - \mu_{r,y,f}^{(\beta)}\right)}{\left(\sigma_{r,y,f}^{(\beta)}\right)^{2}} - 1
```

``` math
a_{r,f}^{(\beta)} = \mu_{r,y,f}^{(\beta)} \cdot \kappa_{r,f}^{(\beta)}
```

``` math
b_{r,f}^{(\beta)} = \left(1 - \mu_{r,y,f}^{(\beta)}\right) \cdot \kappa_{r,f}^{(\beta)}
```

``` math
P\left(\tilde{\beta}_{r,y,f}\right) = \frac{\Gamma\left(a_{r,f}^{(\beta)} + b_{r,f}^{(\beta)}\right)}{\Gamma\left(a_{r,f}^{(\beta)}\right)\Gamma\left(b_{r,f}^{(\beta)}\right)}\tilde{\beta}_{r,y,f}^{a_{r,f}^{(\beta)} - 1}\left(1 - \tilde{\beta}_{r,y,f}\right)^{b_{r,f}^{(\beta)} - 1}
```

where
$`\tilde{\beta}_{r,y,f} = \epsilon + (1 - 2\epsilon) \cdot \text{logistic}\left(\beta_{r,y,f}^{}\right)`$
is the numerically stabilized reporting rate with $`\beta_{r,y,f}^{}`$
the logit-scale parameter and $`\epsilon`$ a small constant
(e.g. $`10^{-4}`$). Here $`\sigma_{r,y,f}^{(\beta)}`$ is the literal
standard deviation of the Beta distribution and must satisfy
$`\sigma_{r,y,f}^{(\beta)} < \sqrt{\mu_{r,y,f}^{(\beta)}\left(1 - \mu_{r,y,f}^{(\beta)}\right)}`$
to ensure $`\kappa_{r,f}^{(\beta)} > 0`$. Note that when
$`\sigma_{r,y,f}^{(\beta)}`$ is large relative to
$`\mu_{r,y,f}^{(\beta)}(1 - \mu_{r,y,f}^{(\beta)})`$, $`a`$ and $`b`$
approach zero and the Beta density becomes U-shaped, placing mass near 0
and 1. In this case, numerical instability can occur during
optimization, which is why $`\tilde{\beta}_{r,y,f}`$ is squished away
from the boundaries via the $`\epsilon`$ transformation.

#### Process Error Penalties

In addition to priors, penalties are also utilized to aid in the
estimation of process errors (either penalized likelihood or integrating
random effects via Laplace Approximation are possible). Currently,
process errors can be specified to arise for initial age deviations,
recruitment, fishing mortality, discard mortality rate, fishery and
survey selectivity, and movement.

##### Initial Age Deviations

To estimate non-equilibrium initial age deviations, multiplicative
deviations can be specified:

``` math
\ell\left( \epsilon_{p,r,i}^{\text{Init}} \right) = \frac{1}{\sqrt{2\pi\sigma_{\text{Init}}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{p,r,i}^{\text{Init}} - \mu^{\text{Init}} \right)^{2}}{2\sigma_{\text{Init}}^{2}} \right)
```

where deviations arise from a normal distribution with variance
$`\sigma_{\text{Init}}^{2}`$. The penalty mean $`\mu^{\text{Init}}`$ is
set by `InitDevs_pen_center`: under `"fixed"` (the default) it is the
asserted prior mean (zero, or the bias-corrected
$`-\sigma_{\text{Rec}}^{2}/2`$ under the bias ramp), constraining both
the level and the spread of the deviations. Under `"own_mean"` it is the
mean of the estimated deviations themselves, so only their spread is
penalized and their level is left free; that is what a sum of squares
about the series’ own mean amounts to, and it matches the convention of
assessments whose deviation vectors are constrained to sum to zero. The
penalty carries its own weight $`\lambda^{\text{InitRec}}`$
(`Wt_Init_Rec`), separate from the recruitment deviation weight, since
the two penalties are dimensioned differently.

When the deviations are sex-specific (`InitDevs_sex_spec = "est_all"`),
each sex’s curve is penalized, while sexes sharing one curve are
penalized once (the shared parameter is not invoked). The `"own_mean"`
center is then pooled over every penalized cell across ages and sexes,
so the sexes share a single estimated level.

A deviation shared across cells of the parameter array through the map,
across regions (`InitDevs_spec = "est_shared_r"`,
`RecDevs_spec = "est_shared_pop_r"`), across sexes, or across the ages
past the observed range (`init_age_devs_shared`), is one parameter and
carries one penalty. Each cell the penalty visits takes $`1/n`$ of it,
with $`n`$ the number of such cells holding that level, so a series
shared by three regions is not penalized three times. Under the bias
ramp the center of the initial age deviation on age $`a`$ is the ramp
read at that age’s own birth year,
$`\mu_{a}^{\text{Init}} = -\sigma_{\text{Rec}}^{2}/2 \cdot b_{1 - a}`$,
the deviation index $`1 - a`$ lying before the first model year, which
is how a ramp defined on calendar years treats the years before the
model starts.

Sex-specific curves can additionally be compared to one another
(`Use_init_sex_pen`), a separate statement about how far apart the
sexes’ initial age structures may sit rather than about how variable
each is:

``` math
\ell\left( \epsilon_{p,r,i,s}^{\text{Init}} - \epsilon_{p,r,i,1}^{\text{Init}} \right) = \frac{1}{\sqrt{2\pi\sigma_{\text{InitSex}}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{p,r,i,s}^{\text{Init}} - \epsilon_{p,r,i,1}^{\text{Init}} \right)^{2}}{2\sigma_{\text{InitSex}}^{2}} \right), \quad s > 1
```

over the same ages $`i`$ the initial-age penalty covers, entering the
objective unweighted with $`\sigma_{\text{InitSex}}`$ held as data
rather than a parameter.

##### Recruitment Deviations

Annual recruitment deviations can also be specified, where
multiplicative deviations are assumed:

``` math
\ell\left( \epsilon_{p,r,y}^{\text{Rec}} \right) = \frac{1}{\sqrt{2\pi}\sigma_{\text{Rec}}^{2}\,}\exp\left( - \frac{\left( \epsilon_{p,r,y}^{\text{Rec}} - \mu^{\text{Rec}} \right)^{2}}{2\sigma_{\text{Rec}}^{2}} \right)
```

with variance $`\sigma_{\text{Rec}}^{2}`$. As for the initial age
deviations, `RecDevs_pen_center` selects the penalty mean
$`\mu^{\text{Rec}}`$: `"fixed"` (default) centers on the asserted prior
mean (zero, or the bias-corrected
$`-\sigma_{\text{Rec}}^{2}/2 \cdot b_y`$ when the bias ramp is active),
whereas `"own_mean"` centers on the mean of the estimated deviations,
penalizing only their spread. Under `"own_mean"` the level of the
deviations is unpenalized and must be pinned elsewhere (an $`R_0`$
prior, or a fixed deviation), or the likelihood is flat along it; it
also cannot be combined with the bias ramp, whose $`-\sigma^2/2`$ offset
is meaningless once the mean is estimated rather than asserted.

Three further refinements act on which deviations are penalized and how
strongly. First, only estimated deviations are penalized: cells mapped
off by hand (via `map$ln_RecDevs`) are excluded from the penalty as well
as from estimation. Second, the recruitment weight
$`\lambda^{\text{Rec}}`$ (`Wt_Rec`) may be a per-deviation array rather
than a scalar, so individual deviations can be down-weighted or removed
from the penalty (weight zero) while remaining estimated. This is how a
stock-recruit relationship is fit over a chosen window of years while
recruitment stays effectively free elsewhere, and is distinct from
`dont_est_recdev_last`, which removes the deviations themselves so
recruitment reverts to the deterministic prediction.

Third, an optional penalty on the level of the recruitment series itself
(`Use_rec_level_pen`) can be applied, separately from the deviation
penalty:

``` math
\ell\left( \log R_{p,r,y} \right) = \frac{1}{\sqrt{2\pi\sigma_{\text{RecLevel}}^{2}}\,}\exp\left( - \frac{\left( \log R_{p,r,y} - \mu^{\text{RecLevel}} \right)^{2}}{2\sigma_{\text{RecLevel}}^{2}} \right)
```

evaluated over a chosen set of years (`rec_level_pen_yrs`), with
$`\mu^{\text{RecLevel}}`$ either the mean of the log series
(`"own_mean"`, penalizing only its variability) or zero (`"fixed"`).
Under a stock-recruit relationship the deviations are residuals about
the predicted curve, so a model that also wants to keep the realized
recruitment series from wandering has nowhere else to say so; this
penalty is that second, independent statement, and reproduces the
recruitment regularity penalties several existing assessments carry.

##### Initial Recruitment Offset

When the initial age structure is built from its own recruitment level
(`use_rinit = 1`, so $`\mu_{p}^{\text{RecInit}}`$ is separate from
$`\mu_{p}^{\text{Rec}}`$), that level can be penalized toward the
recruitment level rather than left free (`Use_rinit_pen`):

``` math
L_{\text{RecInit}} = - \sum_{p = 1}^{n_{p}}\log\left\lbrack \phi\left( \log\mu_{p}^{\text{RecInit}} - \log\mu_{p}^{\text{Rec}};\ 0,\ \sigma_{\text{RecInit}} \right) \right\rbrack
```

The two levels are otherwise only linked through the data, and the
initial age structure is often the thinnest part of it, so the penalty
is a statement about how far the stock’s equilibrium recruitment may
have sat from the level the modeled years show. An equilibrium
recruitment stands for an average over several years rather than one, so
$`\sigma_{\text{RecInit}}`$ is normally narrower than $`\sigma_{R}`$. A
reasonable choice is $`\sigma_{R}/\left( 1/M - 0.5 \right)`$, the
recruitment variability shrunk by the number of year classes the
equilibrium averages over.

##### Stock-Recruit Residual

When recruitment arises about a mean and a stock-recruit curve is
evaluated alongside the dynamics without generating them (`sr_penalty`;
see the Stock-Recruit Curve as a Penalty section above), the curve
enters the objective function through the log residual
$`\xi_{p,r,y}^{\text{SR}} = \log R_{p,r,y} - \log{\widehat{R}}_{p,r,y}`$
between the realized recruitment and the curve’s prediction:

``` math
\ell\left( \xi_{p,r,y}^{\text{SR}} \right) = \frac{1}{\sqrt{2\pi\sigma_{\text{SR}}^{2}}\,}\exp\left( - \frac{\left( \xi_{p,r,y}^{\text{SR}} \right)^{2}}{2\sigma_{\text{SR}}^{2}} \right)
```

The residual is centered on zero, and $`\sigma_{\text{SR}}`$
(`sr_pen_sigma`) is a fixed input rather than an estimated parameter, so
the contribution is a sum of squares weighted by
$`1/\left( 2\sigma_{\text{SR}}^{2} \right)`$ up to a constant; a
template carrying a weight $`w`$ on the squared residuals corresponds to
$`\sigma_{\text{SR}} = 1/\sqrt{2w}`$. The contribution is summed over
the years named by `sr_pen_yrs`, which defaults to every year that has a
lagged spawning biomass behind it, that is, all but the first
$`RecLag`$. Those early years take the equilibrium substitution
described in the Recruitment Processes section rather than
$`{SSB}_{p,r,y - RecLag}`$, so their residual would not be a
stock-recruit residual and naming one is an error rather than a silent
fallback. Years outside the window keep their recruitment deviations
estimated and contribute nothing to this penalty, which is how a
restricted stock-recruit window is expressed.

This penalty and the recruitment deviation penalty are separate
statements about the same series. The deviation penalty acts on
$`\epsilon_{p,r,y}^{\text{Rec}}`$, an estimated parameter, about its own
center, whereas this one acts on the difference between two derived
quantities, so the deviations stay free to depart from the curve at a
cost governed by $`\sigma_{\text{SR}}`$. It sits with the process error
penalties rather than the priors above because it constrains a model
quantity rather than asserting knowledge about a parameter.

The penalty is a statement about the objective function alone and is not
a data-generating process. Simulation, whether self-testing or
closed-loop, generates recruitment from the mean parameter and its
deviations exactly as the mean-recruitment equation specifies; no
simulated recruitment is drawn from the curve, and
$`{\widehat{R}}_{p,r,y}`$ plays no part in the simulated population
dynamics.

##### Fishing Mortality Deviations

Fishing mortality deviations assume multiplicative deviations about a
mean rate. One of three process error structures can be specified via
`Fdev_model`: independent (`"iid"`), random walk (`"rw"`), or
first-order autoregressive (`"ar1"`). In all three cases, the penalty is
only evaluated in region $`r`$, season $`\tau`$, and fleet $`f`$
combinations with observed catch (i.e.,
$`\text{UseCatch}_{r,y,\tau,f} = 1`$ or any
$`\text{UseCatch\_pop}_{p,r,y,\tau,f} = 1`$).

###### IID

``` math
\ell\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} \right) = \frac{1}{\sqrt{2\pi\sigma_{r,\tau,f,\text{Fsh}}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} - \mu^{\text{Fsh}} \right)^{2}}{2\sigma_{r,\tau,f,\text{Fsh}}^{2}} \right)
```

Fishing mortality deviations are assumed to arise from a normal
distribution with variance $`\sigma_{r,\tau,f,\text{Fsh}}^{2}`$. The
penalty mean $`\mu^{\text{Fsh}}`$ is set by `Fdev_pen_center` (iid model
only): `"fixed"` (default) centers on zero, constraining both the level
and the spread of the deviations, while `"own_mean"` centers on the mean
of the estimated deviations, penalizing only their spread. Under the
mean-plus-deviations parameterization the level of fishing mortality is
already carried by $`\mu_{r,\tau,f}^{\text{Fsh}}`$, so `"own_mean"`
avoids penalizing it twice; note it also leaves $`\mu^{\text{Fsh}}`$ and
the deviations’ level mutually unidentified unless one of them is fixed,
which `ln_F_mean_spec = "fix"` does (free annual log-F, see the fishing
mortality parameterization above).

###### Random Walk

Catch-active years need not be contiguous under the random walk (a
fishery may close for several years and reopen later). Let $`y'`$ denote
the previous catch-active year for a given region, season, and fleet,
and $`d = y - y'`$ the number of elapsed years between them ($`d = 1`$
when catch is available every year). The first catch-active year is
initialized with a large, diffuse variance; every subsequent
catch-active year follows a random walk about the previous active year’s
value, with variance inflated by the elapsed gap $`d`$:

``` math
\ell\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} \right) = \left\{ \begin{matrix}
\frac{1}{\sqrt{2\pi \cdot 5^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} \right)^{2}}{2 \cdot 5^{2}} \right),\ y\ \text{is the first catch-active year} \\
\frac{1}{\sqrt{2\pi d\sigma_{r,\tau,f,\text{Fsh}}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} - \epsilon_{r,y',\tau,f}^{\text{Fsh}} \right)^{2}}{2d\sigma_{r,\tau,f,\text{Fsh}}^{2}} \right),\ \text{otherwise} \\
\end{matrix} \right.\ 
```

When catch is available every year ($`d = 1`$ throughout), this reduces
exactly to a standard single-step random walk. When years are closed
(e.g., a fishery closure), inflating the variance by $`d`$ gives exactly
the same marginal distribution that would be obtained by estimating
deviations for the closed years and integrating them out without
actually estimating them, so no deviation parameters exist for closed
years.

###### AR1

The AR1 form additionally estimates a correlation parameter,
$`\rho_{r,\tau,f,\text{Fsh}} \in (-1,1)`$ (from an unconstrained
parameter `Fdev_rho`, transformed via $`\rho = 2/(1+e^{-2x}) - 1`$). As
with the random walk, catch-active years need not be contiguous. The
first catch-active year is drawn from the process’s stationary marginal
distribution, and every subsequent catch-active year follows an AR1
transition over the elapsed gap $`d`$ since the previous active year:

``` math
\ell\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} \right) = \left\{ \begin{matrix}
\frac{1}{\sqrt{2\pi\sigma_{r,\tau,f,\text{Fsh}}^{2}/\left(1 - \rho_{r,\tau,f,\text{Fsh}}^{2}\right)}\,}\exp\left( - \frac{\left(1 - \rho_{r,\tau,f,\text{Fsh}}^{2}\right)\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} \right)^{2}}{2\sigma_{r,\tau,f,\text{Fsh}}^{2}} \right),\ y\ \text{is the first catch-active year} \\
\frac{1}{\sqrt{2\pi V_{d}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,\tau,f}^{\text{Fsh}} - \rho_{r,\tau,f,\text{Fsh}}^{d}\epsilon_{r,y',\tau,f}^{\text{Fsh}} \right)^{2}}{2V_{d}} \right),\ \text{otherwise} \\
\end{matrix} \right.\ 
```

where
``` math
V_{d} = \sigma_{r,\tau,f,\text{Fsh}}^{2}\sum_{i = 0}^{d - 1}{\rho_{r,\tau,f,\text{Fsh}}^{2i}} = \sigma_{r,\tau,f,\text{Fsh}}^{2}\frac{1 - \rho_{r,\tau,f,\text{Fsh}}^{2d}}{1 - \rho_{r,\tau,f,\text{Fsh}}^{2}}
```
is the exact variance of the sum of the $`d`$ intervening (unestimated)
innovations that would have occurred during the closed years, and
$`\rho_{r,\tau,f,\text{Fsh}}^{d}`$ is the corresponding decay of the
mean across the same gap. As with the random walk, this reduces exactly
to the standard single-step AR1 transition when $`d = 1`$.

##### Discard Mortality Rate Deviations

Discard mortality rate deviations are penalized analogously on the logit
scale:

``` math
\ell\left( \epsilon_{r,y,\tau,f}^{\delta} \right) = \frac{1}{\sqrt{2\pi{\sigma_{r,\tau,f}^{2}}_{\delta}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,\tau,f}^{\delta} \right)^{2}}{2{\sigma_{r,\tau,f}^{2}}_{\delta}} \right)
```

where $`{\sigma_{r,\tau,f}^{2}}_{\delta}`$ is the variance of the
discard mortality rate deviations. The penalty is only applied in years
and fleets where discard data are available.

##### Fishery and Survey Selectivity

A variety of process error parameterizations can be specified for
fishery and survey selectivity. Across all parameterizations,
multiplicative deviations are assumed. In the most basic case, iid
deviations can be assumed to vary about a parameter on a given
selectivity functional form:

``` math
\ell\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} \right) = \frac{1}{\sqrt{2\pi\sigma_{r,i,s,j,Sel}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} \right)^{2}}{2\sigma_{r,i,s,j,Sel}^{2}} \right)
```

where $`\epsilon_{r,y,i,s,j}^{\text{Sel}}`$ are selectivity deviations
about a given parameter for region $`r`$, year $`y`$, parameter $`i`$,
sex $`s`$, and fleet $`j`$. Deviations are assumed to have a mean of 0
and a variance of $`\sigma_{r,i,s,j,Sel}^{2}`$, constrained by a normal
distribution.

Extending the iid case, random walk selectivity deviations can also be
specified about a given parameter, assuming a normal distribution:

``` math
\ell\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} \right) = \left\{ \begin{matrix}
\frac{1}{\sqrt{2\pi \cdot \sigma_{0,j}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y = 1,i,s,j}^{\text{Sel}} \right)^{2}}{2\sigma_{0,j}^{2}} \right),\ if\ y = 1 \\
\frac{1}{\sqrt{2\pi\sigma_{r,i,s,j,Sel}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} - \epsilon_{r,y - 1,i,s,j}^{\text{Sel}} \right)^{2}}{2\sigma_{r,i,s,j,Sel}^{2}} \right),\ if\ y > 1 \\
\end{matrix} \right.\ 
```

where process error deviations for the first year are initialized with
their own standard deviation $`\sigma_{0,j}`$ (`sel_rw_init_sigma`),
which defaults to a large value ($`5`$) that leaves the first year
effectively unconstrained. Setting it to `NA` instead starts the walk at
zero under the walk’s own $`\sigma_{Sel}`$, which is what a first
difference taken against a selectivity of one amounts to and makes the
first year’s deviation as smooth as every later step. Following the
first year, process error deviations follow a random walk process with a
mean conditional on the previous year’s value
($`\epsilon_{r,y - 1,i,s,j}^{\text{Sel}}`$) and a variance of
$`\sigma_{r,i,s,j,Sel}^{2}`$.

Each fleet’s process error contribution is additionally scaled by a
per-fleet weight $`\lambda_{j}^{\text{SelPE}}`$ (`fishsel_pe_wt`,
`retsel_pe_wt`, `srvsel_pe_wt`; default $`1`$) before entering the joint
likelihood. A weight of zero removes the process error penalty for that
fleet entirely while its deviations remain estimated, which is how a
model reproduces assessments that let selectivity deviations float
subject only to explicit smoothness penalties rather than a
distributional assumption.

The bin-override deviations described in the selectivity section carry
their own process error (`cont_tv_sel_bin_devs`), either iid or a random
walk of exactly the forms above, with their own estimated $`\sigma`$ per
overridden bin and their own first-year standard deviation
(`sel_bin_devs_rw_init_sigma`). An overridden bin that is free of the
functional form is thereby still smoothed in time when asked to be.

In addition to being constrained by a normal distribution, both iid and
random walk cases have an optional additional smoothness penalty
applied:

``` math
P_{SelYrSmooth} = \sum_{r = 1}^{n_{r}}{\sum_{j = 1}^{n_{j}}{\sum_{s = 1}^{n_{s}}{\sum_{y = 2}^{n_{y}}{\sum_{i = 1}^{n_{i}}\left( \epsilon_{r,y,i,s,j}^{\text{Sel}}- \epsilon_{r,y - 1,i,s,j}^{\text{Sel}} \right)^{2}}}}}
```

Additionally, semi-parametric deviations can also be specified. In
total, there are three options that can be utilized, two of which allow
age, year, and cohort correlations, while one allows for only age or
length and year correlations. In the case where age, year, and cohort
correlations are specified (note that this is only possible when
age-based selectivity is specified), marginal stationary variance and a
conditional non-stationary variance can be invoked. The following
equations describe the conditional variance version:

``` math
\mathbf{\epsilon}_{r,s,j}^{\text{Sel}}\mathbf{=}vec\left( \epsilon_{r,y,a,s,j}^{\text{Sel}} \right)
```

where we vectorize the selectivity deviations across its year and age
dimensions. These deviations are then assumed to arise from a
multivariate normal distribution (or Gaussian Markov Random Field) with
a covariance matrix ($`\mathbf{\Sigma} = \mathbf{Q}^{- 1}`$) determined
by:

``` math
\begin{matrix}
\text{Q} = \left( \text{I} - \left( \text{B} \right)^{T} \right)\mathbf{\Omega}\left( \text{I} - \text{B} \right) \\
\text{diag}\left( \mathbf{\Omega} \right) = \sigma_{r,s,j,Sel}^{- 2} \\
\end{matrix}
```

Here, $`\text{I}`$ is an identity matrix and $`\mathbf{\Omega}`$ is a
diagonal matrix that determines the variance of the multivariate normal
process. $`\text{B}`$ is a square matrix representing the partial effect
of $`\mathbf{\epsilon}_{r,s,j}^{\text{Sel}}`$ on preceding ages and/or
years, governed by partial correlation coefficients for ages, years, and
cohorts. To demonstrate the formulation of $`\text{B}`$, a simplified
example is provided with rows representing ages
$`a\ \epsilon\ \{ 1,2\}`$ and columns representing years
$`y\ \epsilon\ \{ 1,2\}`$. In this example, $`\text{B}`$ is a
$`4\ \times\ 4`$ matrix, where both the rows and columns represent
combinations of age and year:

``` math
\begin{matrix}
\mathbf{B} = \begin{bmatrix}
1 & 0 & 0 & 0 \\
\rho_{y} & 0 & 0 & 0 \\
\rho_{a} & 0 & 0 & 0 \\
\rho_{c} & \rho_{a} & \rho_{y} & 0 \\
\end{bmatrix} \\
\end{matrix}
```

where $`\rho_{y}`$, $`\rho_{a}`$, and $`\rho_{c}`$ are parameters
describing the partial autocorrelation among years within a given age,
among ages within a given year, and years within a cohort, respectively.
The multivariate likelihood is then defined as:

``` math
\ell\left( \mathbf{\epsilon}_{r,s,j}^{\text{Sel}} \right) = \frac{\left| \mathbf{Q} \right|^{1/2}}{(2\pi)^{\frac{n}{2}}}\exp\left( - \frac{1}{2}\left( \mathbf{\epsilon}_{r,s,j}^{\text{Sel}} \right)^{T}\mathbf{Q}\left( \mathbf{\epsilon}_{r,s,j}^{\text{Sel}} \right) \right)
```

If age or length and year correlations are specified (i.e., a
two-dimensional autoregressive structure), a multivariate normal
likelihood is similarly assumed, but the covariance structure of this
process is defined as:

``` math
\text{Q}^{- 1} = \frac{\sigma_{r,s,j,Sel}^{2}}{\left( 1 - \rho_{y} \right)^{2}\left( 1 - \rho_{b} \right)^{2}}\text{R}_{y} \otimes \text{R}_{b}
```

where $`\rho_{y}`$ and $`\rho_{b}`$ are correlation coefficients across
years and bins, respectively.

##### Selectivity Smoothness Penalties

A set of six penalty terms, evaluated directly on a fleet’s realized
selectivity-at-bin-at-year surface rather than on any particular
selectivity parameterization, can be independently weighted and applied
to any selectivity functional form.

The dome-shape penalty discourages the selectivity curve from decreasing
across adjacent bins within a year (i.e., encourages flat-topped or
asymptotic rather than dome-shaped curves, when desired), applied only
where an actual decrease occurs:

``` math
P_{SelSmoothDome} = \sum_{r = 1}^{n_{r}}{\sum_{j = 1}^{n_{j}}{\sum_{s = 1}^{n_{s}}{\sum_{y}{\sum_{b}{\left\lbrack \max\left( \log\left( {Sel}_{r,y,b,s,j} \right) - \log\left( {Sel}_{r,y,b + 1,s,j} \right),\ 0 \right) \right\rbrack^{2}}}}}}
```

The bin (age or length) curvature penalty is a second-difference
smoothness penalty across bins, normalized by the number of fitted bins
$`n_{b}^{\text{fit}}`$:

``` math
P_{SelSmoothBinCurve} = \frac{1}{n_{b}^{\text{fit}}}\sum_{r = 1}^{n_{r}}{\sum_{j = 1}^{n_{j}}{\sum_{s = 1}^{n_{s}}{\sum_{y}{\sum_{b}\left( \log\left( {Sel}_{r,y,b + 1,s,j} \right) - 2\log\left( {Sel}_{r,y,b,s,j} \right) + \log\left( {Sel}_{r,y,b - 1,s,j} \right) \right)^{2}}}}}
```

A related, unconditional first-difference penalty across bins where both
increases and decreases contribute, unlike the dome-shape penalty above
which is normalized the same way:

``` math
P_{SelSmoothBinDiff} = \frac{1}{n_{b}^{\text{fit}}}\sum_{r = 1}^{n_{r}}{\sum_{j = 1}^{n_{j}}{\sum_{s = 1}^{n_{s}}{\sum_{y}{\sum_{b}\left( \log\left( {Sel}_{r,y,b,s,j} \right) - \log\left( {Sel}_{r,y,b + 1,s,j} \right) \right)^{2}}}}}
```

Inter-annual variation is penalized with a first-difference penalty
across years, normalized by the number of fitted years
$`n_{y}^{\text{fit}}`$:

``` math
P_{SelSmoothYrDiff} = \frac{1}{n_{y}^{\text{fit}}}\sum_{r = 1}^{n_{r}}{\sum_{j = 1}^{n_{j}}{\sum_{s = 1}^{n_{s}}{\sum_{b}{\sum_{y}\left( \log\left( {Sel}_{r,y,b,s,j} \right) - \log\left( {Sel}_{r,y - 1,b,s,j} \right) \right)^{2}}}}}
```

and inter-annual smoothness with an analogous second-difference penalty
across years:

``` math
P_{SelSmoothYrCurve} = \frac{1}{n_{y}^{\text{fit}}}\sum_{r = 1}^{n_{r}}{\sum_{j = 1}^{n_{j}}{\sum_{s = 1}^{n_{s}}{\sum_{b}{\sum_{y}\left( \log\left( {Sel}_{r,y + 1,b,s,j} \right) - 2\log\left( {Sel}_{r,y,b,s,j} \right) + \log\left( {Sel}_{r,y - 1,b,s,j} \right) \right)^{2}}}}}
```

Finally, because some selectivity forms (e.g. the bicubic spline) have
no built-in scale identifiability constraint (a uniform per-year shift
in log-selectivity trades off exactly against that year’s fishing
mortality), a mean-centering penalty regularizes the per-year mean of
log-selectivity toward zero:

``` math
P_{SelSmoothMeanCenter} = \sum_{r = 1}^{n_{r}}{\sum_{j = 1}^{n_{j}}{\sum_{s = 1}^{n_{s}}{\sum_{y}\left\lbrack \frac{1}{n_{b}^{\text{fit}}}\sum_{b}{\log\left( {Sel}_{r,y,b,s,j} \right)} \right\rbrack^{2}}}}
```

Each of the six terms above ($`P_{SelSmoothDome}`$,
$`P_{SelSmoothBinCurve}`$, $`P_{SelSmoothBinDiff}`$,
$`P_{SelSmoothYrDiff}`$, $`P_{SelSmoothYrCurve}`$,
$`P_{SelSmoothMeanCenter}`$) is scaled by its own
independently-specified weight before being added to the joint negative
log-likelihood, allowing each to be turned on or off and tuned
separately. In code, these six weights use a `smooth_` prefix
(e.g. `smooth_bin_curve`, `smooth_yr_diff`) rather than referencing the
bicubic spline specifically, since, as described above, they apply to
any selectivity form.

Several generalizations of these penalties are available, specified per
fleet (a single specification is shared by every fleet, or an unnamed
list gives each fleet its own):

- Per-year weights. Every weight $`w`$ may be a vector with one value
  per model year, $`w_{y}`$, so a penalty can act only in the years
  where selectivity is allowed to change, or act with a different
  strength in each year. Years with $`w_{y} = 0`$ are skipped entirely.
  With `normalize = FALSE` (below), a per-year first-difference weight
  of $`w_{y} = 1/(2\sigma_{y}^{2})`$ makes $`P_{SelSmoothYrDiff}`$
  exactly the negative log-kernel of a random walk with year-specific
  standard deviation $`\sigma_{y}`$, expressed as a penalty on the
  realized curve rather than on deviation parameters.
- Bin range. `bin_range` restricts the bins the penalties act over to
  $`\lbrack b_{lo}, b_{hi} \rbrack`$, either one range shared by every
  term or a named list giving each term its own. This is how a shape
  penalty (e.g., the dome or curvature penalty) is confined to the older
  ages where the curve is expected to flatten, without constraining the
  ascending limb.
- Normalization. `normalize` controls whether the bin-wise weights are
  divided by the number of penalized bins and the year-wise weights by
  the number of years (the default, `TRUE`); it can likewise be set per
  term. Turn it off when weights are calibrated as explicit variances,
  as in the random walk correspondence above.
- First-year reference. The year first-difference walk has no
  predecessor in its first penalized year, so that year is normally left
  unpenalized. `yr_diff_ref` supplies one: the first penalized year is
  held toward a reference log-selectivity vector, anchoring an otherwise
  free series to a known selectivity before the data begin.

##### Growth Deviations

Growth carries two deviation surfaces, and both are scored by the same
process error machinery the selectivity deviations use, so the
vocabulary above carries over unchanged.

A time-varying growth parameter has one deviation per year, a surface
one column wide. Under `"iid"` each is drawn independently and under
`"rw"` they form a random walk:

``` math
\begin{matrix}
\delta_{p,r,y,k,s} \sim N\left( 0, \sigma_{k}^{2} \right) & \texttt{iid} \\[1ex]
\delta_{p,r,y,k,s} \sim N\left( \delta_{p,r,y - 1,k,s}, \sigma_{k}^{2} \right) & \texttt{rw}
\end{matrix}
```

with one $`\sigma_{k}`$ per varying parameter $`k`$, held in the first
stream of `growth_pe_pars` and estimated under
`growth_tv_sigma_spec = "est"`, and the first year of a walk given its
own standard deviation (`growth_rw_init_sigma`). Only the years named in
`growth_tv_years` carry a deviation; the rest are held at zero, so a
parameter can be constant early in a series and varying once the data
can support it.

The semi-parametric surface $`\varepsilon_{p,r,y,a,s}`$ runs over years
and ages, so all five structures are available to it: `"iid"` and `"rw"`
as above, `"2dar1"` for a separable first-order autoregression over ages
and years, and `"3dmarg"`/`"3dcond"` for the three-dimensional Gaussian
Markov random field over age, year and cohort described for selectivity,
on the marginal or conditional variance. The correlated forms are the
ones that make a growth surface estimable in practice: the deviations
are not identified year by year and age by age from length data alone,
and it is the correlation that lets neighboring ages and years share
information. Its hyperparameters live in the second stream of
`growth_pe_pars`, in the same slots the selectivity forms use, and the
surface can be restricted to the ages and years the length data inform
through `growth_semipar_ages` and `growth_semipar_years`.

Both enter the objective unweighted, reported as $`L_{\text{GrowthTV}}`$
and $`L_{\text{GrowthSemipar}}`$.

##### Movement

Time-varying movement is introduced through process error deviations,
$`\epsilon^{\text{Move}}`$, which modify baseline movement parameters.
The interpretation of these deviations depends on the movement
formulation, but their stochastic structure is shared.

###### General Structure

Movement deviations are assumed to be independent and normally
distributed:

``` math
\epsilon^{\text{Move}}_{p,r,r',y,\tau,a,s}
\sim
N\left(0, \sigma^2_{p,r,\tau,a,s,\text{Move}}\right)
```

where $`\sigma_{p,r,\tau,a,s,\text{Move}}`$ may be shared across
dimensions depending on the selected process error model.

Only valid origin-destination pairs (i.e., adjacent regions) are
assigned deviations.

###### Unstructured Markov Movement

For multinomial logit movement, deviations enter additively in logit
space:

``` math
\omega_{p,r,k,y,\tau,a,s}
=
\omega_{p,r,k,\tau,s}
+
\epsilon^{\text{Move}}_{p,r,k,y,\tau,a,s}
```

Thus, time variation is expressed as year-specific perturbations around
a mean logit, and movement probabilities are obtained via the softmax
transform.

###### CTMC Movement

For CTMC movement, deviations act on the transition rates rather than
logits. Specifically, deviations are applied multiplicatively to the
off-diagonal diffusion terms:

``` math
D_{p,r \to r',y,\tau,a,s}
=
\bar{D}_{p,r \to r',y^,\tau,a,s}
\cdot
\exp\left(
\epsilon^{\text{Move}}_{p,r,r',y,\tau,a,s}
\right)
\quad \text{for } r \ne r'
```

where: - $`\bar{D}_{p,r \to r',y^,\tau,a,s}`$ is the baseline diffusion
rate (constructed from covariates and parameters, with year lookups
capped at $`y^ = \min(y, n_{\text{yrs}})`$), -
$`\epsilon^{\text{Move}}_{p,r,r',y,\tau,a,s}`$ is the deviation applied
on the log scale.

This formulation implies that: - deviations are log-multiplicative on
movement rates, - $`\exp(\epsilon^{\text{Move}})`$ acts as a
proportional scaling factor, - time variation persists into projection
years even when baseline covariates are held fixed.

Deviations are applied only to off-diagonal elements (i.e., actual
transitions), and the diagonal of the generator matrix is recomputed to
preserve mass balance.

###### Likelihood for Deviations

The movement process error contribution to the log-likelihood can be
written explicitly as:

``` math
\ell_{\text{Move}}
=
\sum_{p,r,r',y,\tau,a,s}
\left[
-\frac{1}{2}\log\left(2\pi \sigma^2_{p,r,\tau,a,s,\text{Move}}\right)
-
\frac{
\left(\epsilon^{\text{Move}}_{p,r,r',y,\tau,a,s}\right)^2
}{
2\sigma^2_{p,r,\tau,a,s,\text{Move}}
}
\right]
```

where the summation is taken over all valid origin-destination pairs
(i.e., $`r \neq r'`$ and adjacency$`(r,r') = 1`$), and over all indices
of population ($`p`$), year ($`y`$), season ($`\tau`$), age ($`a`$), and
sex ($`s`$).

###### Variance Structures

Different process error models specify how $`\sigma`$ is shared across
dimensions. These correspond to IID assumptions over subsets of:

- population ($`p`$),
- year ($`y`$),
- season ($`\tau`$),
- age ($`a`$),
- sex ($`s`$).

For example: - IID across years: $`\sigma_{p,r}`$, - IID across years
and ages: $`\sigma_{p,r,a}`$, - Fully stratified:
$`\sigma_{p,r,\tau,a,s}`$.

These structures control the degree of temporal and demographic
heterogeneity in movement variability.

#### Joint Likelihood

Lastly, the joint likelihood to be minimized represents the sum of all
observational likelihood components, priors, and penalties defined
above:

``` math
\begin{matrix}
\text{Joint Likelihood} = \sum\text{Observation Likelihoods} + \sum\text{Priors} + \sum\text{Penalties} \\
\end{matrix}
```

Note that some of these components may be zero (i.e., if no priors are
used) depending on the configuration of the model. The recruitment and
initial age deviation penalties enter with their own separate weights
($`\lambda^{\text{Rec}}`$, possibly per deviation, and
$`\lambda^{\text{InitRec}}`$), while the recruitment level penalty and
the stock-recruit residual penalty
$`\sum_{p,r,y}^{}{-\log\ell\left( \xi_{p,r,y}^{\text{SR}} \right)}`$,
when enabled, each enter as their own additional component carrying no
weight of their own.

## References

Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H., Bell, B., 2016. TMB:
Automatic Differentiation and Laplace Approximation. J. Stat. Soft. 70.
<https://doi.org/10.18637/jss.v070.i05>

Mace, P.M., Doonan, I.J., 1988. A Generalized Bioeconomic Simulation
Model for Fish Population Dynamics. MAFFish, N.Z. Ministry of
Agriculture and Fisheries.

McGarvey, R., Feenstra, J.E., 2002. Estimating rates of fish movement
from tag recoveries: conditioning by recapture. Can. J. Fish. Aquat.
Sci. 59, 1054-1064. <https://doi.org/10.1139/f02-080>

Methot, R.D., Taylor, I.G., 2011. Adjusting for bias due to variability
of estimated recruitments in fishery assessment models. Can. J. Fish.
Aquat. Sci. 68, 1744-1760. <https://doi.org/10.1139/f2011-092>

Monnahan, C.C., 2024. Toward good practices for Bayesian data-rich
fisheries stock assessments using a modern statistical workflow.
Fisheries Research 275, 107024.
<https://doi.org/10.1016/j.fishres.2024.107024>

Thorson, J.T., Johnson, K.F., Methot, R.D., Taylor, I.G., 2017.
Model-based estimates of effective sample size in stock assessment
models using the Dirichlet-multinomial distribution. Fisheries Research
192, 84-93. <https://doi.org/10.1016/j.fishres.2016.06.005>
