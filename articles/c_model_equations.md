# Description of Model Equations

The Stochastic Population over Regional Components (`SPoRC`) model is a
generalized integrated population model written in `R-TMB` (R bindings
for Template Model Builder; Kristensen et al., 2016) that supports age,
sex, and spatially-structured dynamics. Currently, `SPoRC` assumes that
population dynamics operate following an annual time step, and in
sequential order:

1.  Recruitment and tag releases initially occur (tag releases only
    occur if tagging data are used),
2.  Markovian movement of individuals then follow (movement only occurs
    in the spatial model), and
3.  Total mortality, chronic tag loss (if applicable), and ageing
    processes occur.

These processes are modeled across four primary model partitions: region
($r$), year ($y$), age ($a$ and $a_{+}$, where $a_{+}$ is the plus
group), and sex ($s$). In single-region and/or single-sex models, these
equations generally reduce by setting $r = 1$ and/or $s = 1$. In
general, the same equations are used for both simulation and estimation.

### Process Equations

#### Population Initialization

In `SPoRC`, three primary methods exist to initialize the equilibrium
population of the model. The first method derives the equilibrium
population using the following process:

$$\begin{matrix}
{N_{r,a,s}^{\prime} = \mu^{\text{Rec}}\exp\left( - (a - 1) \cdot \left( Z_{r,a,s}^{\prime} \right) \right)\psi_{r,y = 1,s}\zeta_{r},\quad{\text{for}\mspace{6mu}}2 \leq a < a_{+}} \\

\end{matrix}$$

$$Z_{r,a,s}^{\prime} = \text{Natmort}_{r,y = 1,a,s} + \mu_{r,f = 1}^{\text{Fsh}}\text{Fmort}^{\text{InitProp}}\text{Sel}_{r,y = 1,a,s,f = 1}^{\text{Fsh}}$$

where:

- $N_{r,a,s}^{\prime}$ are the equilibrium numbers-at-age,
- $\mu^{\text{Rec}}$ is a global recruitment parameter (either virgin or
  mean recruitment depending on the parameterization of stock
  recruitment dynamics),
- $Z_{r,a,s}^{\prime}$ is the initial instantaneous total mortality
  rate,
- $\text{Natmort}_{r,y = 1,a,s}$ is the instantaneous natural mortality
  rate,
- $\mu_{r,f = 1}^{\text{Fsh}}$ is the mean fishing mortality rate from
  the first fleet ($f$ denotes fishery fleet),
- $\text{Fmort}^{\text{InitProp}}$ is a user-defined value describing
  the proportion of the mean fishing mortality from the first fleet
  applied during the initialization stage,
- $\text{Sel}_{r,y = 1,a,s,f = 1}^{\text{Fsh}}$ is the fishery
  selectivity-at-age for the first fleet,
- $\psi_{r,y,s}$ describes the recruitment sex-ratio,
- $\zeta_{r}$ apportions the global recruitment parameter across regions
  (estimated using a multinomial logit transform to ensure proportions
  sum to one).

The plus group ($a_{+}$) of the initial population is then computed as:

$$\begin{matrix}
{N_{r,a_{+},s}^{\prime} = N_{r,a_{+} - 1,s}^{\prime}\frac{\exp\left( - Z_{r,a = a_{+} - 1,s}^{\prime} \right)}{1 - \exp\left( Z_{r,a = a_{+},s}^{\prime} \right)}} \\

\end{matrix}$$

However, this scalar geometric series solution assumes that the plus
group accumulates in a closed system. Therefore, when movement dynamics
are present, this solution does not correctly accumulate individuals
into the plus group.

To address this, additional methods are provided to explicitly
incorporate movement dynamics into the plus group calculation. In
particular, the initial population can be derived by iterating the
population to equilibrium. Here, the total instantaneous mortality in
each region is defined similarly as in the first method. An exponential
decay model is used to first initialize the age structure at the first
iteration:

$$\begin{matrix}
{N_{r,a,s}^{\prime} = \left\{ \begin{matrix}
{\mu^{\text{Rec}}\,\psi_{r,y = 1,s}\,\zeta_{r},} & {{\text{if}\mspace{6mu}}a = 1} \\
{\mu^{\text{Rec}}\,\psi_{r,y = 1,s}\,\zeta_{r}\,\exp\left( - \sum\limits_{a = 1}^{n_{a}}Z_{r,a,s}^{\prime} \right),} & {{\text{if}\mspace{6mu}}a > 1} \\
 & 
\end{matrix} \right.\ } \\

\end{matrix}$$

The initialized age structure is then iterated forward to equilibrium by
applying recruitment, movement, and mortality and ageing processes in
order (see Population Projection section for equations).

While the iterative method correctly accumulates the plus group when
movement is present, it can be computationally inefficient due to the
potentially large number of iterations required to propagate the initial
population. Therefore, SPoRC enables users to compute the plus group
using the matrix formulation of the geometric series, which correctly
accounts for movement processes. Individuals are first initialized
through recruitment processes, after which movement, mortality, and
ageing processes are applied sequentially (up to the number of ages in
the model). The population (or cohort) is then projected forward to the
penultimate age ($N_{r,a_{+} - 1,s}^{\prime}$), which is used to derive
the plus group. Specifically, the penultimate age is projected forward
once more:

$$\mathbf{X}_{s} = \left( \left( \mathbf{N}_{a_{+} - 1,s}^{\prime} \right)^{T}\mathbf{M}_{y = 1,a_{+} - 1,s} \right)\text{diag}\left( \exp\left( - \mathbf{Z}_{a = a_{+} - 1,s}^{\prime} \right) \right)$$

where $\mathbf{M}$ is a first-order Markov matrix representing movement,
and $\mathbf{X}_{s}$ represents the culmination of processes applied to
the penultimate age. A transition matrix ${(\mathbf{G}}_{s}$) is the
constructed to represent the combined effects of survival and movement
on the plus group:

$$\mathbf{G}_{s} = \text{diag}\left( \exp\left( - \mathbf{Z}_{a = a_{+},s}^{\prime} \right) \right)\left( \mathbf{M}_{y = 1,a_{+},s} \right)^{T}$$

The plus group solution incorporating movement is then given by:

$$\mathbf{N}_{a_{+},s}^{\prime} = \left( {\mathbf{I} -}\mathbf{G}_{s} \right)^{- 1}\mathbf{X}_{s}$$

When only a single region is modeled or no movement occurs (i.e., an
identity matrix), the matrix formulation simplifies to the standard
scalar geometric series solution.

Following the definition of equilibrium age structure, initial age
deviations can be applied:

$$\begin{matrix}
{N_{r,y = 1,a \neq 1,s} = N\prime_{r,a \neq 1,s}\text{exp}\left( \epsilon_{r,i}^{\text{Init}} \right)} \\

\end{matrix}$$

where $N_{r,y = 1,a \neq 1,s}$ are represents the numbers-at-age in the
first model year except recruits ($a \neq 1$). These values can then be
treated as a stochastic process by applying multiplicative lognormal
deviations $\epsilon_{r,i}^{\text{Init}}$ to the initial equilibrium age
structure. Note that the index $i$ is introduced because users can
determine whether initial age deviations are estimated up to the
penultimate age class, or across all classes including the plus group.

#### Recruitment Processes

In the current iteration of `SPoRC`, two stock recruitment
parameterizations can be specified. Recruitment can be specified to
arise about a mean parameter ($\mu^{\text{Rec}}$):

$$\begin{matrix}
{N_{r,y,a = 1,s} = \mu^{\text{Rec}}\exp\left( \epsilon_{r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)\psi_{r,y,s}\zeta_{r}} \\

\end{matrix}$$

where $\epsilon_{r,y}$ are annual, lognormally distributed recruitment
deviations with a lognormal bias correction term
($\frac{\sigma_{\text{Rec}}^{2}}{2}b_{y}$), with $b_{y}$ representing
the bias correction ramp from Methot and Taylor (2011). Recruitment can
also be specified to arise from a Beverton-Holt stock recruitment
function to invoke density-dependent population dynamics, which follows
the steepness parameterization (Mace and Doonan, 1988). Localized
density-dependent recruitment is defined as:

$$\begin{matrix}
{N_{r,y,a = 1,s} = \frac{4\mu^{\text{Rec}}\zeta_{r}h_{r}{SSB}_{r,y - RecLag}}{\left( 1 - h_{r} \right)SSB0_{r} + 5\left( h_{r} - 1 \right){SSB}_{r,y - RecLag}}\exp\left( \epsilon_{r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)\psi_{r,y,s}} \\

\end{matrix}$$

while global density-dependent recruitment can be defined as:

$$\begin{matrix}
{N_{r,y,a = 1,s} = \frac{4\mu^{\text{Rec}}\zeta_{r}h\sum\limits_{r}^{}{SSB}_{r,y - RecLag}}{(1 - h)\sum\limits_{r}^{}{SSB0_{r}} + 5(h - 1)\sum\limits_{r}^{}{SSB}_{r,y - RecLag}}\exp\left( \epsilon_{r,y}^{\text{Rec}} - \frac{\sigma_{\text{Rec}}^{2}}{2}b_{y} \right)\psi_{r,y,s}} \\

\end{matrix}$$

where $\mu^{\text{Rec}}$ under this parameterization is the virgin
unfished recruitment, $h_{r}$ (or $h$) is the steepness parameter
representing the fraction of $\mu^{\text{Rec}}\zeta_{r}$ that would be
produced when at 20% of $SSB0_{r}$ (or $\sum_{r}^{}{SSB0_{r}}$). The
steepness parameter is constrained to be between values of 0.2 and 1 and
are estimated in bounded logit space. $SSB0_{r}$ is a derived variable
that represents the unfished spawning stock biomass (i.e., derived by
multiplying the unfished SSB per recruit by
$\mu^{\text{Rec}}\zeta_{r}$). ${SSB}_{r,y - RecLag}$ is the spawning
stock biomass, which is the product of numbers-at-age, spawning
weight-at-age, and maturity-at-age for females:

$$\begin{matrix}
{SSB_{r,y} = \sum\limits_{a = 1}^{a_{+}}{N_{r,y,a,s = 1}W_{r,y,a,s = 1}^{spawn}{Mat}_{r,y,a,s = 1}}} \\

\end{matrix}$$

Note that $RecLag$ denotes the delay between spawning and when recruits
enter the population. Thus, if $RecLag > y$, `SPoRC` utilizes $SSB0_{r}$
instead of ${SSB}_{r,y - RecLag}$ to compute deterministic recruitment
($SSB0_{r}$ is further decremented if
$\mu_{r,f = 1}^{\text{Fsh}}\text{Fmort}^{\text{InitProp}} \neq 0$).

#### Population Projection

Following recruitment processes, the population can then be projected
forward. In the context of the spatial model, Markovian movement
dynamics are first applied:

$$\begin{matrix}
{\mathbf{N}_{y,a,s} = \left( \mathbf{N}_{y,a,s} \right)^{T}\mathbf{M}_{y,a,s}} \\

\end{matrix}$$

Here, $\mathbf{M}_{y,a,s}$ is a first-order Markov matrix representing
movement. In a single-region case, no movement is applied (i.e.,
$\mathbf{M}_{y,a,s}$ is an implied identity matrix). For each year, age,
and sex combination, the movement matrix specifies bulk-transfer
coefficients, with parameters transformed through a multinomial logit to
ensure that proportions within each origin region sum to one. Following
movement processes, the population undergoes mortality and ageing:

$$\begin{matrix}
{N_{r,y + 1,a + 1,s} = N_{r,y,a,s}\exp\left( - Z_{r,y,a,s} \right),\quad{\text{for}\mspace{6mu}}1 \leq a < a_{+}} \\

\end{matrix}$$

$$N_{r,y + 1,a,s} = N_{r,y,a - 1,s}\exp\left( - Z_{r,y,a - 1,s} \right) + N_{r,y,a,s}\exp\left( - Z_{r,y,a,s} \right),\quad{\text{for}\mspace{6mu}}a = a_{+}$$

$Z_{r,y,a,s}$ denotes the annual total instantaneous mortality rate and
is defined as the combination of natural mortality
($\text{Natmort}_{r,y,a,s}$) and fishing mortality
($\text{Fmort}_{r,y,f}$):

$$\begin{matrix}
{Z_{r,y,a,s} = \text{Natmort}_{r,y,a,s} + \sum\limits_{f}^{}\text{Sel}_{r,y,a,s,f}^{\text{Fsh}}\text{Fmort}_{r,y,f}} \\

\end{matrix}$$

The annual instantaneous fishing mortality rate is then defined as:

$$\begin{matrix}
{\text{Fmort}_{r,y,f} = \mu_{r,f}^{\text{Fsh}}\text{exp}\left( \epsilon_{r,y,f}^{\text{Fsh}} \right)} \\

\end{matrix}$$

where $\text{Fmort}_{r,y,f}$ is parameterized based on lognormal
deviations ($\epsilon_{r,y,f}^{\text{Fsh}})$ about a mean fishing
mortality parameter for a given fishery fleet
($\mu_{r,f}^{\text{Fsh}}$).

#### Movement Processes

Movement processes can be parameterized as either an unstructured Markov
process (discrete-time Markov) or as a Continuous-time Markov chain
(CTMC) process. Movement parameterized as an unstructured Markov process
is estimated using a multinomial logit link function, with
$n_{r} \times n_{r} - 1$ free parameters:

$$M_{r,k,y,a,s} = \frac{\exp\left( \omega_{r,k,y,a,s} \right)}{\sum\limits_{k}\exp\left( \omega_{r,k,y,a,s} \right)}$$

where the reference region $k = 1$ is set as
$\omega_{r,k = 1,y,a,s} = 0$. Under this parameterization, movement
fractions can be estimated independently for each stratum $(y,a,s)$, or
grouped into blocks to reduce the number of parameters. While parameter
blocking can help reduce the parameter load, estimating all pairwise
movement probabilities under the unstructured Markov approach can still
result in a large number of parameters.

Alternatively, movement can be specified as a continuous-time Markov
chain (CTMC) process, which decomposes into diffusive and taxis
components. Diffusive processes represent undirected movement of
individuals, while taxis processes represent directed movement toward
more preferred habitat. This CTMC movement parameterization is governed
by an adjacency matrix ($A_{r,k}$), which defines neighboring regions
that can receive individuals within a given time step. Diffusive
processes are given by:

$${\dot{D}}_{r,k,y,a,s} = \begin{cases}
{\frac{e^{2\theta}}{V_{r}},} & {{\text{if}\mspace{6mu}}A_{r,k} = 1{\mspace{6mu}\text{and}\mspace{6mu}}r \neq k,} \\
{- \sum\limits_{j \neq r}{\dot{D}}_{r,j,y,a,s},} & {{\text{if}\mspace{6mu}}r = k,} \\
{0,} & \text{otherwise.}
\end{cases}$$

where ${\dot{D}}_{r,k,y,a,s}$ represents diffusion, $\theta$ is the log
diffusion rate, $V_{r}$ scales the log diffusion rate, such that smaller
regions have higher diffusion rates, and $j$ in the second equation
indexes all destinations except the source to ensure that the rows of
the matrix sum to 0, thereby conserving abundance. Taxis processes
(preference) can then be written as:

$${\dot{P}}_{r,k,y,a,s} = \begin{cases}
{h_{k,y,a,s} - h_{r,y,a,s},} & {{\text{if}\mspace{6mu}}A_{r,k} = 1{\mspace{6mu}\text{and}\mspace{6mu}}r \neq k,} \\
{- \sum\limits_{j \neq r}{\dot{P}}_{r,j,y,a,s},} & {{\text{if}\mspace{6mu}}r = k,} \\
{0,} & \text{otherwise.}
\end{cases}$$ Here, ${\dot{P}}_{r,k,y,a,s}$ represents the taxis
(preference) component of movement. Equation 1 determines local
differences in the habitat preference function, $\mathbf{h}_{y,a,s}$,
while equation 2 ensures that the rows of the matrix sum to 0,
conserving abundance.

Habitat preference can be defined flexibly as a combination of linear
effects and basis splines:

$$h_{r,y,a,s} = \sum\limits_{k = 1}^{n_{k}}\beta_{r,k}W_{r,k,y,a,s}$$

where $\beta_{r,k}$ are the estimated effects (incorporating linear or
spline effects), and $W_{r,k,y,a,s}$ is the design matrix.

Diffusive and taxis processes can then be combined to construct a
generator matrix:

$${\dot{Q}}_{r,k,y,a,s} = {\dot{D}}_{r,k,y,a,s} + {\dot{P}}_{r,k,y,a,s}$$

Here, ${\dot{Q}}_{r,k,y,a,s}$ represents instantaneous movement rates.
The generator matrix ${\dot{\mathbf{Q}}}_{y,a,s}$ is Metzler, such that
all off-diagonal elements satisfy:

$${\dot{Q}}_{r,k,y,a,s} \geq 0\quad{\text{for}\mspace{6mu}}r \neq k.$$

Given that movement in SPoRC is defined sequentially, the instantaneous
movement matrix can then be converted to annual movement fractions using
the matrix exponential:

$$\mathbf{M}_{y,a,s} = \exp({\dot{\mathbf{Q}}}_{y,a,s}\,\Delta t)$$

where $\Delta t$ is the duration of the movement interval.

### Observation Equations

#### Fishery Observation Model

The fishery observation model describes the expected catch-at-age,
catch-at-length, catch (in units of biomass), and fishery indices.
Expected catch-at-age ($C_{r,y,a,s,f}^{a}$) for a given fishery fleet is
calculated using Baranov’s catch equation:

$$\begin{matrix}
{C_{r,y,a,s,f}^{a} = \frac{\text{Fmort}_{r,y,f}\text{Sel}_{r,y,a,s,f}^{\text{Fsh}}}{Z_{r,y,a,s}}N_{r,y,a,s}\left\lbrack 1 - \exp\left( - Z_{r,y,a,s} \right) \right\rbrack} \\

\end{matrix}$$

To track length-based dynamics, catch-at-length ($C_{r,y,l,s,f}^{l}$)
can then be derived using the following equation:

$$\begin{matrix}
{C_{r,y,l,s,f}^{l} = \mathbf{A}_{r,y,s}^{l}{\mathbf{C}^{\mathbf{a}}}_{r,y,s,f}} \\

\end{matrix}$$

where $\mathbf{A}_{r,y,s}^{l}$ is a user-defined size-age transition
matrix to convert ages to lengths. Expected catch
($\text{Catch}_{r,y,f}$) can then be computed as either abundance or
biomass. For catch in units of abundance, this is given by the summation
of the expected catch-at-age across ages and sexes:

$$\begin{matrix}
{\text{Catch}_{r,y,f} = \sum\limits_{a}^{a_{+}}{\sum\limits_{s}^{n_{s}}C_{r,y,a,s,f}^{a}}} \\

\end{matrix}$$

For catch in biomass, this is given by the summation of expected
catch-at-age, multiplied by their respective fishery weights-at-age
($W_{r,y,a,s,f}^{fish}$):

$$\begin{matrix}
{\text{Catch}_{r,y,f} = \sum\limits_{a}^{a_{+}}{\sum\limits_{s}^{n_{s}}C_{r,y,a,s,f}^{a}}W_{r,y,a,s,f}^{fish}} \\

\end{matrix}$$

Similarly expected fishery indices ($\text{FshIdx}_{r,y,f}$) can be
computed as either abundance-based or biomass-based. For abundance-based
fishery indices, this is computed as:

$$\begin{matrix}
\begin{matrix}
{\text{FshIdx}_{r,y,f} = q_{r,y,f}^{\text{Fsh}}\sum\limits_{a}^{a_{+}}{\sum\limits_{s}^{n_{s}}N_{r,y,a,s}}\exp\left( - Z_{r,y,a,s} \right)\text{Sel}_{r,y,a,s,f}^{\text{Fsh}}} \\

\end{matrix} \\

\end{matrix}$$

where $q_{r,y,f}^{\text{Fsh}}$ is the catchability coefficient for a
given fishery fleet. Biomass-based fishery indices are computed as:

$$\begin{matrix}
\begin{matrix}
{\text{FshIdx}_{r,y,f} = q_{r,y,f}^{\text{Fsh}}\sum\limits_{a}^{a_{+}}{\sum\limits_{s}^{n_{s}}N_{r,y,a,s}}\exp\left( - Z_{r,y,a,s} \right)\text{Sel}_{r,y,a,s,f}^{\text{Fsh}}W_{r,y,a,s,f}^{fish}} \\

\end{matrix} \\

\end{matrix}$$

#### Survey Observation Model

Likewise, the survey observation model describes the expected survey
catch-at-age, survey catch-at-length, and survey indices. Expected
survey catch-at-age ($I_{r,y,a,s,sf}^{a}$) is calculated as follows:

$$\begin{matrix}
{I_{r,y,a,s,sf}^{a} = N_{r,y,a,s}\exp\left( - Z_{r,y,a,s}t^{srv} \right)\text{Sel}_{r,y,a,s,sf}^{\text{Srv}}} \\

\end{matrix}$$

where subscript $sf$ denotes a given survey fleet and
$\text{Sel}_{r,y,a,s,sf}^{\text{Srv}}$ is the survey selectivity-at-age
pattern. Expected survey catch-at-length ($I_{r,y,l,s,sf}^{l}$) is given
by:

$$\begin{matrix}
{I_{r,y,l,s,sf}^{l} = \mathbf{A}_{r,y,s}^{l}{\mathbf{I}^{\mathbf{a}}}_{r,y,s,f}} \\

\end{matrix}$$

Survey indices ($\text{SrvIdx}_{r,y,sf}$) can be computed as either
abundance-based or biomass-based. Abundance-based survey indices are
calculated as:

$$\begin{matrix}
{\text{SrvIdx}_{r,y,sf} = q_{r,y,sf}^{\text{Srv}}\sum\limits_{a}^{a_{+}}{\sum\limits_{s}^{n_{s}}I_{r,y,a,s,sf}^{a}}} \\

\end{matrix}$$

while biomass-based indices are computed as:

$$\begin{matrix}
{\text{SrvIdx}_{r,y,sf} = q_{r,y,sf}^{\text{Srv}}\sum\limits_{a}^{a_{+}}{\sum\limits_{s}^{n_{s}}I_{r,y,a,s,sf}^{a}}W_{r,y,a,s,sf}^{srv}} \\

\end{matrix}$$

Here, $W_{r,y,a,s,sf}^{srv}$ is the weight-at-age for a given survey,
$q_{r,y,sf}^{\text{Srv}}$ represents the survey catchability coefficient
in the aforementioned equations.

For survey catchability, environmental linkage can be specified:

$$q_{r,y,sf}^{\text{Srv}} = q_{r,sf}^{\text{Srv}}\exp\left( \mathbf{x}^{T}{{\mathbf{β}} +}\sum\limits_{m}^{}{\iota_{m}p_{m}\left( z_{r,y,sf} \right)} \right)$$

where $q_{r,sf}^{\text{Srv}}$ is the base survey catchability (i.e.,
intercept), $\mathbf{x}$ is a matrix of covariates, $\mathbf{β}$ is a
vector of regression coefficients, $\iota_{m}p_{m}$ are orthogonal
polynomial coefficients along with its basis functions, and $z_{r,y,sf}$
are the covariates for which a polynomial term is assumed for.

### Tagging Observation Model

The tagging observation model tracks tag cohorts ($T_{r,y,a,s}^{k}$) by
the combination of release region and release year ($k$) and follows a
Brownie tag attrition framework. Tag cohorts are tracked for a
pre-defined maximum duration (maximum tag liberty; $n_{L}$), after which
calculations for the tag cohort are no longer computed. This is done to
aid with computation time. Additionally, tag cohorts also have the
option to undergo a fraction of mortality in the year of release, such
that if a given cohort is released during the middle of the year,
mortality processes are discounted. If mortality discounting is
specified, movement does not occur in the year of release. In general,
the process dynamics for the tagged cohort mimics those specified for
the overall population. Immediately following release, tag cohorts are
decremented by an initial tag-induced mortality rate:

$$\begin{matrix}
{T_{r,y,a,s}^{k} = T_{r,y,a,s}^{k}\exp( - \tau)} \\

\end{matrix}$$

where $\tau$ is the initial tag-induced mortality rate. If tagged
cohorts are released at the beginning of a calendar year, Markovian
movement occurs (movement does not occur otherwise):

$$\begin{matrix}
{\mathbf{T}_{y,a,s}^{k} = \left( \mathbf{T}_{y,a,s}^{k} \right)^{\text{T}}\mathbf{M}_{y,a,s}} \\

\end{matrix}$$

Mortality and ageing of the tagged cohort then occurs:

$$\begin{matrix}
\begin{matrix}
{T_{r,y + 1,a + 1,s}^{k} = T_{r,y,a,s}^{k}\exp\left( - Z_{r,y,a,s}^{\text{Tag}} \right),\quad{\text{for}\mspace{6mu}}2 \leq a < a_{+}} \\
{T_{r,y + 1,a,s}^{k} = T_{r,y,a - 1,s}^{k}\exp\left( - Z_{r,y,a - 1,s}^{\text{Tag}} \right) + T_{r,y,a,s}^{k}\exp\left( - Z_{r,y,a,s}^{\text{Tag}} \right),\quad{\text{for}\mspace{6mu}}a = a_{+}} \\

\end{matrix} \\

\end{matrix}$$

where tagged cohorts follow an exponential mortality model, with
accumulation of individuals in the plus-group. Total mortality for the
tagged cohort ($Z_{r,y,a,s}^{\text{Tag}}$) is specified as:

$$\begin{matrix}
{Z_{r,y,a,s}^{\text{Tag}} = \left( \kappa + \text{NatMort}_{r,y,a,s} + \sum\limits_{f}^{}\text{Sel}_{r,y,a,s,f}^{\text{Fsh}}\text{Fmort}_{r,y,f} \right)} \\

\end{matrix}$$

where $\kappa$ is a parameter describing chronic tag loss (i.e., annual
tag shedding). Similar to computations for catch-at-age, tag recaptures
are calculated using a modified version of Baranov’s catch equation:

$$\begin{matrix}
{\text{Recap}_{r,y,a,s}^{k} = \beta_{r,y}\frac{\sum\limits_{f}^{}\text{Sel}_{r,y,a,s,f}^{\text{Fsh}}}{Z_{r,y,a,s}^{\text{Tag}}}T_{r,y,a,s}^{k}\left\lbrack 1 - \exp\left( - Z_{r,y,a,s}^{\text{Tag}} \right) \right\rbrack} \\

\end{matrix}$$

$\beta_{r,y}$ represents a tag reporting rate parameter that can vary
across space and time, which is estimated in logit space such that it is
constrained between $\lbrack 0,1\rbrack$.

### Fishery and Survey Selectivity

In the following descriptions of selectivity, we omit subscripts for
sexes and fleets for brevity, although note that the equations remain
specific to those model partitions. Several approaches are available for
parameterizing fishery and survey selectivity. Selectivity can be
defined as either age- or length-based. For age-based selectivity, an
age vector is applied directly with a chosen functional form. For
length-based selectivity, a length vector is used to compute
selectivity, which is then combined with a size–age transition matrix
through a dot product to obtain selectivity-at-age.

$$\begin{matrix}
{{\mathbf{S}\mathbf{e}}\mathbf{l}^{a} = \left( \mathbf{A}_{r,y,s}^{l} \right)^{\text{T}}{\mathbf{S}\mathbf{e}}\mathbf{l}^{l}} \\

\end{matrix}$$

Given the age-based nature of the model, selectivity-at-age is utilized
for all subsequent calculations. Various functional forms can also be
specified for selectivity processes. In the following we use the
subscript $b$ to denote a generalized bin number.

Two forms of logistic selectivity can be specified. The first form is
defined as:

$$\begin{matrix}
{{Sel}_{b} = \frac{1}{1 + \exp\left\lbrack - k\left( b - b^{50} \right) \right\rbrack}} \\

\end{matrix}$$

where $k$ determines the slope/steepness of the logistic curve and
$b^{50}$ is the bin-at-50% selection. The second form can be expressed
as:

$$\begin{matrix}
{{Sel}_{b} = \frac{1}{1 + 19^{(\frac{b^{50} - b}{b^{95}})}}} \\

\end{matrix}$$

Here, $b^{50}$ is also the bin-at-50% selection and $b^{95}$ is the
bin-at-95% selection. Beyond the specification of flat-topped
selectivity, `SPoRC` also allows for dome-shaped selectivity. In
particular, a reparametrized gamma function can be specified:

$$\begin{matrix}
{p = 0.5\left\lbrack \sqrt{b^{\max} + 4\gamma^{2}} - b^{\max} \right\rbrack} \\

\end{matrix}$$

$${Sel}_{b} = \left( \frac{b}{b^{\max}} \right)^{\frac{b^{\max}}{p}}\exp\left( \frac{b^{\max} - b}{p} \right)$$

In this parameterization, $p$ is a derived power parameter, $\gamma$ is
the shape parameter that describes the steepness of the descending limb,
and $b^{\max}$ describes the bin-at-maximum selection. Dome-shaped
selectivity can also be expressed as a power function:

$$\begin{matrix}
{{Sel}_{b} = \frac{1}{b^{\phi}}} \\

\end{matrix}$$

with $\phi$ being a power parameter that determines the descending limb
of the curve (larger values are steeper). The last dome-shaped
selectivity form that can be specified includes a 6 parameter (denoted
as ${\widehat{p}}_{1}$ through ${\widehat{p}}_{6}$) double normal
functional form with the following transformations applied:

$$\begin{matrix}
{p_{1} = \ \min(b) + \left\lbrack \max(b) - min(b) \right\rbrack \cdot \frac{1}{1 + \exp\left( {\widehat{p}}_{1} \right)}} \\

\end{matrix}$$

$$p_{2} = \ p_{1} + 1 + \frac{0.99 + max(b) - p_{1} - 1}{1 + exp\left( {\widehat{p}}_{2} \right)}$$

$$p_{3} = \exp\left( {\widehat{p}}_{3} \right),\ \ p_{4} = \exp\left( {\widehat{p}}_{4} \right)$$

$$p_{5} = \frac{1}{1 + \exp\left( {\widehat{p}}_{5} \right)},\ \ p_{6} = \frac{1}{1 + \exp\left( {\widehat{p}}_{6} \right)}$$

$p_{1}$ represents the beginning of the plateau, $p_{2}$ describes the
width of the plateau, $p_{3}$ and $p_{4}$ are parameters controlling the
ascending and descending widths, and $p_{5}$ and $p_{6}$ control the
selectivity at the first and last bins. The double normal function is
then defined with a series of functions:

$$\begin{matrix}
{{asc}_{b} = \exp\left( - \frac{\left( b - p_{1} \right)^{2}}{p_{3}} \right)} \\
{{ascscaled}_{b} = p_{5} + \left( 1 - p_{5} \right) \cdot {asc}_{b}} \\
{{desc}_{b} = \exp\left( - \frac{\left( b - p_{2} \right)^{2}}{p_{4}} \right)} \\
{stj = exp\left( - \frac{\left( 40 - p_{2} \right)^{2}}{p_{4}} \right)} \\
{{descscaled}_{b} = 1 + \left( p_{6} - 1 \right) \cdot \frac{{desc}_{b} - 1}{stj - 1}} \\
{join1_{b} = \frac{1}{1 + \exp\left( - 20 \cdot \frac{b - p_{1}}{1 + \left| b - p_{1} \right|} \right)}} \\
{join2_{b} = \frac{1}{1 + \exp\left( - 20 \cdot \frac{b - p_{2}}{1 + \left| b - p_{2} \right|} \right)}} \\

\end{matrix}$$

and are joined together as:

$$\begin{matrix}
{{Sel}_{b} = {asc}_{b} \cdot \left( 1 - join1_{b} \right) + join1_{b} \cdot \left\lbrack \left( 1 - join2_{b} \right) + {descscaled}_{b} \cdot join2_{b} \right\rbrack} \\

\end{matrix}$$

The double normal functional form is incredibly flexible and is able to
reduce to both flat-topped and dome-shaped selectivity forms, depending
on the values of the parameters.

In addition to the functional forms that can be specified to describe
selectivity processes, several options exist to specify continuous
time-varying processes. In particular, options to specify time-varying
parametric selectivity and time-varying semi-parametric selectivity are
available. To illustrate, if logistic selectivity is specified and
parametric deviations are invoked, the following expression is used:

$$\begin{matrix}
{{Sel}_{y,b} = \frac{1}{1 + \exp\left\lbrack - k_{y}\left( b - b_{y}^{50} \right) \right\rbrack}} \\

\end{matrix}{k_{y} = k \cdot \exp\left( \epsilon_{y,1}^{Sel} \right)}{b_{y}^{50} = b^{50} \cdot exp\left( \epsilon_{y,2}^{Sel} \right)}$$

where the parameters of the logistic form are allowed to vary over time.
In the context of semi-parametric selectivity, the following equation is
used:

$$\begin{matrix}
{{Sel}_{y,b}^{\prime} = \frac{1}{1 + \exp\left( - k\left\lbrack b - b^{50} \right\rbrack \right)}\exp\left( {\epsilon_{y,b}}^{Sel} \right)} \\

\end{matrix}{{Sel}_{y,b} = \frac{{Sel}_{y,b}^{\prime}}{mean\left( {\mathbf{S}\mathbf{e}}\mathbf{l}^{\mathbf{\prime}} \right)}}$$

where deviations are placed about the parametric form and selectivity
values are mean standardized to aid with interpretability. Importantly,
note that if age-based selectivity is specified and continuous
deviations are estimated, these deviations are placed on the ages
themselves. In the case where length-based selectivity is specified,
these deviations are placed on the length bins. Further details on how
selectivity deviations arise can be found in the “Selectivity Process
Error” section of this document.

## Likelihoods

Currently, `SPoRC` incorporates data likelihood components for the
following data sources:

1.  fishery catches,
2.  fishery indices,
3.  fishery age compositions,
4.  fishery length compositions,
5.  survey indices,
6.  survey age compositions,
7.  survey length compositions, and
8.  tagging data.

The total likelihood (objective function) is the sum of the individual
likelihood contributions from these data sources along with priors and
penalties, where the objective function is minimized using a non-linear
optimization algorithm to estimate model parameters.

### Observation Likelihoods

#### Fishery Catches

Fishery catches can be fit using a lognormal likelihood. The
log-likelihood for observed catch,
$L\left( \log\left( \text{ObsCatch}_{r,y,f} \right) \right)$ is defined
as:

$$\begin{matrix}
{L\left( \log\left( \text{ObsCatch}_{r,y,f} \right) \right) =} \\

\end{matrix}$$

$$\lambda_{\text{ObsCatch}_{r,y,f}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsCatch}_{r,,f}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsCatch}_{r,y,f} \right) - log\left( \text{Catch}_{r,y,f} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsCatch}_{r,y,f}}^{2}} \right)$$

Here, $\lambda_{\text{ObsCatch}_{r,y,f}}$ is the likelihood weight,
$\text{ObsCatch}_{r,y,f}$ is the observed catch, $\text{Catch}_{r,y,f}$
is the predicted catch, and $\sigma_{\text{ObsCatch}_{r,y,f}}^{2}$ is
the variance of catch on the log scale.

#### Fishery and Survey Indices

Fishery indices can also be fit assuming a lognormal likelihood. The
log-likelihood for observed fishery indices is:

$$\begin{matrix}
{L\left( \log\left( \text{ObsFshIdx}_{r,y,f} \right) \right) =} \\

\end{matrix}$$

$$\lambda_{\text{ObsFshIdx}_{r,y,f}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsFshIdx}_{r,y,f}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsFshIdx}_{r,y,f} \right) - log\left( \text{FshIdx}_{r,y,f} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsFshIdx}_{r,y,f}}^{2}} \right)$$

where $\lambda_{\text{ObsFshIdx}_{r,y,f}}$ controls the weight of
fishery indices to the objective function, $\text{ObsFshIdx}_{r,y,f}$
represents the observed fishery indices, and
$\sigma_{\text{ObsFshIdx}_{r,y,f}}^{2}$ denotes the variance of the
fishery index.

Likewise, survey indices can be fit assuming a lognormal likelihood
given by:

$$\begin{matrix}
{L\left( \log\left( \text{ObsSrvIdx}_{r,y,sf} \right) \right) =} \\

\end{matrix}$$

$$\lambda_{\text{ObsSrvIdx}_{r,y,sf}} \cdot \frac{1}{\sqrt{2\pi\sigma_{\text{ObsSrvIdx}_{r,y,sf}}^{2}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{ObsSrvIdx}_{r,y,sf} \right) - log\left( \text{SrvIdx}_{r,y,sf} \right) \right\rbrack^{2}}{2\sigma_{\text{ObsSrvIdx}_{r,y,sf}}^{2}} \right)$$

$\lambda_{\text{ObsSrvIdx}_{r,y,sf}}$ is the likelihood weight applied
to survey indices, $\text{ObsSrvIdx}_{r,y,sf}$ are the observed survey
indices, and $\sigma_{\text{ObsSrvIdx}_{r,y,sf}}^{2}$ indicates the
variance of the survey index.

#### Fishery and Survey Compositions

Several options for fitting composition data are available in `SPoRC`.
These include the multinomial, the Dirichlet-multinomial, and the
logistic-normal likelihoods. In the case of the multinomial likelihood,
the following expression is used:

$$\begin{matrix}
{L\left( \text{ObsCompositionData}_{r,y,j} \right) =} \\

\end{matrix}{\lambda_{\text{ObsCompositionData}_{r,y,j}} \cdot \text{ISS}_{r,y,j} \cdot \prod\limits_{b = 1}^{n_{b}}E_{r,y,b,j}^{O_{r,y,b,j}}}$$

where subscript $j$ is used to indicate a fishery or survey fleet and
the $b$ subscript generically indicates a bin number.
$\lambda_{\text{ObsCompositionData}_{r,y,j}}$ are likelihood weights
applied to composition data (can be derived using iterative weighting
methods), $ISS_{r,y,j}$ is the input sample size, $E_{r,y,b,j}$ denotes
the expected composition proportions, and $O_{r,y,b,j}$ are the observed
composition proportions.

If a Dirichlet-multinomial likelihood is assumed, the following
parameterization (linear) is used:

$$\begin{matrix}
{L\left( \text{ObsCompositionData}_{r,y,j} \right) =} \\
{\lambda_{\text{ObsCompositionData}_{r,y,j}} \cdot \frac{\Gamma\left( \text{ISS}_{r,y,j} + 1 \right)}{\text{ISS}_{r,y,j} \cdot O_{r,y,b,j} + 1} \cdot \frac{\Gamma\left( \text{ISS}_{r,y,j}\theta_{r,j} \right)}{\text{Γ(}\text{ISS}_{r,y,j}\text{+}\text{ISS}_{r,y,j}\theta_{r,j}\text{)}} \cdot} \\
{\prod\limits_{b = 1}^{n_{b}}\frac{\Gamma\left( \text{ISS}_{r,y,j}O_{r,y,b,j} + \text{ISS}_{r,y,j}\theta_{r,j}E_{r,y,b,j} \right)}{\text{ISS}_{r,y,j}\theta_{r,j}E_{r,y,b,j}}} \\

\end{matrix}$$

Here, $\theta_{r,j}$ is the overdispersion parameter of the
Dirichlet-multinomial that adjusts the input sample size. The effective
sample size ($\text{ESS}_{r,y,j})$ can then be derived as:

$$\begin{matrix}
{\text{ESS}_{r,y,j} = \frac{1}{1 + \theta_{r,j}\ } + \text{ISS}_{r,y,j}\frac{\theta_{r,j}}{1 + \theta_{r,j}\ }} \\

\end{matrix}$$

Therefore, if $\theta_{r,j}$ was large such that
$\theta_{r,j} \gg \text{ISS}_{r,y,j}$, then
$\left. \text{ESS}_{r,y,j}\rightarrow\text{ISS}_{r,y,j} \right.$. By
contrast, if $\theta_{r,j}$ was estimated at a small value such that
$\theta_{r,j} \ll \text{ISS}_{r,y,j}$, then $\theta_{r,j}$ approximates
the ratio of $\text{ESS}_{r,y,j}$ and $\text{ISS}_{r,y,j}$ (Thorson et
al., 2017).

A multivariate logistic-normal likelihood can also be assumed, which is
given by:

$$\begin{matrix}
{L\left( \text{ObsCompositionData}_{r,y,j} \right) =} \\

\end{matrix}{\frac{1}{(2\pi)^{\frac{B - 1}{2}}|\mathbf{\Sigma}|^{\frac{1}{2}}}\exp\left( - \frac{1}{2}\left\{ \mathbf{O}_{r,y,j}^{\mathbf{\prime}} - \mathbf{E}_{r,y,j}^{\mathbf{\prime}} \right\}^{T}\mathbf{\Sigma}^{- 1}\left\{ \mathbf{O}_{r,y,j}^{\mathbf{\prime}} - \mathbf{E}_{r,y,j}^{\mathbf{\prime}} \right\} \right)}$$

Both $\mathbf{O}_{r,y,j}^{\mathbf{\prime}}$ and
$\mathbf{E}_{r,y,j}^{\mathbf{\prime}}$ are $(B - 1)$ dimensional
vectors, while $\mathbf{\Sigma}$ is a $(B - 1) \times (B - 1)$
covariance matrix (see below for further details).
$\mathbf{O}_{r,y,j}^{\mathbf{\prime}}$ and
$\mathbf{E}_{r,y,j}^{\mathbf{\prime}}$ are derived via an additive
logistic function:

$$\begin{matrix}
{O_{r,y,b,j}^{\prime} = \log\left( O_{r,y, - B,j} \right) - log\left( O_{r,y,B,j} \right)} \\
{E_{r,y,b,j}^{\prime} = \log\left( E_{r,y, - B,j} \right) - log\left( E_{r,y,B,j} \right)} \\

\end{matrix}$$

where $O_{r,y,b,j}^{\prime}$ and $E_{r,y,b,j}^{\prime}$ are transformed
proportions using the last bin $B$ as the reference category. Because
the logarithm of zero is undefined, all untransformed proportions must
be strictly positive. If any observed proportion is zero, both the
observed and corresponding expected values are removed, and the
remaining proportions are renormalized to ensure that they sum to one
before applying the transformation. The covariance matrix of the
logistic-normal likelihood can be specified in various ways. In the
simplest case, the covariance matrix can be assumed to be independent
and identically distributed (iid):

$$\begin{matrix}
{{\mathbf{\Sigma} =}\left( \mathbf{I}_{B} \cdot \theta^{2} \right)_{- B}} \\

\end{matrix}$$

where $\mathbf{I}_{B}$ is a $B \times B$ identity matrix and
$\theta^{2}$ is an estimated overdispersion parameter representing the
variance. The simple iid case can be further extended to incorporate a
one-dimensional lag-1 autoregressive structure:

$$\begin{matrix}
{{\mathbf{\Sigma} =}\left( \mathbf{R}_{B} \cdot \frac{\theta^{2}}{1 - \rho_{B}^{2}} \right)_{- B}} \\
{\left( \mathbf{R}_{B} \right)_{i,j} = \rho_{B}^{|i - j|},\ \ i,j = \ 1,\ \cdots,\ B} \\

\end{matrix}$$

Here, $\mathbf{R}_{B}$ is a $B \times B$ correlation matrix with a lag-1
autoregressive structure, where $\rho_{B}^{|i - j|}$ defines the
correlation across bins. Lastly, if the model is specified to be
sex-structured and sex-composition data are utilized, a two-dimensional
autoregressive structure can be specified:

$$\begin{matrix}
{{\mathbf{\Sigma} =}\left( {\mathbf{R}_{S}{\ \bigotimes\ \mathbf{R}}}_{C} \cdot \frac{\theta^{2}}{\left( 1 - \rho_{S}^{2} \right)\left( 1 - \rho_{C}^{2} \right)} \right)_{- B}} \\

\end{matrix}$$

$$\left( \mathbf{R}_{C} \right)_{i,j} = \rho_{C}^{|i - j|},\ \ i,j = \ 1,\ \cdots,\ C$$

$$\left( \mathbf{R}_{S} \right)_{i,j} = \left\{ \begin{matrix}
{1,\ \ ifi = j} \\
{\rho_{s},\ \ ifi \neq j} \\

\end{matrix} \right.\ ,\ \ i,j = 1,\ldots,n_{s}$$

$\mathbf{R}_{S}$ is a constant correlation matrix dimensioned by
$n_{s} \times n_{s}$for sexes, with off-diagonal elements $\rho_{s}$
controlling the correlation of age/length categories across sexes, while
$\mathbf{R}_{C}$ is a $n_{c} \times n_{c}$ lag-1 autoregressive
correlation structure, where $\rho_{C}^{|i - j|}$ defines the
correlation across age/length categories. $\bigotimes$ denotes the
Kronecker product. This structure implies that similar age or size
categories across sexes are correlated, while ensuring that very
different age or size categories remain uncorrelated.

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

$$\begin{matrix}
{E_{y,b}^{\prime} = \frac{\sum\limits_{r = 1}^{n_{r}}{\sum\limits_{s = 1}^{n_{s}}E_{r,y,b,s}^{\prime}}}{n_{s} \cdot n_{r}}} \\

\end{matrix}{\mathbf{E}_{y} = \mathbf{E}_{y}^{\mathbf{\prime}}\mathbf{\Theta}_{y}}$$

where compositions are summed across regions and sexes and normalized to
sum to one. Ageing error ($\mathbf{\Theta}_{y}$) can then be applied
using standard matrix multiplication. Expected compositions that are
specified as ‘Split’ by sexes and regions are computed as:

$$\begin{matrix}
{E_{y,b,s}^{\prime} = \frac{E_{r,y,b,s}^{\prime}}{\sum\limits_{b = 1}^{n_{B}}E_{r,y,b,s}^{\prime}}} \\

\end{matrix}$$

$$\mathbf{E}_{y,s} = \mathbf{E}_{y,s}^{\mathbf{\prime}}\mathbf{\Theta}_{y}$$

Here, expected compositions sum to one within a given region and sex
combination and ageing error is similarly applied via matrix
multiplication. In the case where expected compositions are specified as
‘Joint’, they are calculated as:

$$\begin{matrix}
{E_{y,b,s}^{\prime} = \frac{E_{r,y,b,s}^{\prime}}{\sum\limits_{s = 1}^{n_{s}}{\sum\limits_{b = 1}^{n_{B}}E_{r,y,b,s}^{\prime}}}} \\

\end{matrix}$$

$$\mathbf{E}_{y} = \left( \mathbf{E}_{y} \right)^{T}\left( \mathbf{I}_{s}{\ \bigotimes\ \mathbf{\Theta}}_{y} \right)$$

where the expected compositions sum to one jointly across bins and
sexes, thus preserving implicit sex-ratio information. Ageing error is
then applied by taking the Kronecker product of a $n_{s}\ x\ n_{s}$
identity matrix with the ageing error matrix, followed by matrix
multiplication.

#### Tagging

`SPoRC` currently allows for various tagging likelihoods, ranging from
the Poisson, Negative Binomial, multinomial, and Dirichlet-multinomial
likelihood. Additionally, `SPoRC` also allows for both release- and
recapture-conditioned dynamics (McGarvey and Feenstra, 2002). The
Poisson tag likelihood is given by:

$$\begin{matrix}
{L\left( {ObsRecap}_{r,y,a,s}^{k} \right) = \lambda_{\text{Tagging}}\frac{\exp\left( - {Recap}_{r,y,a,s}^{k} \right)\left( {Recap}_{r,y,a,s}^{k} \right)^{{ObsRecap}_{r,y,a,s}^{k}}}{{ObsRecap}_{r,y,a,s}^{k}!}} \\

\end{matrix}$$

where ${ObsRecap}_{r,y,a,s}^{k}$ are the observed tag recaptures and
$\lambda_{\text{Tagging}}$ is the likelihood weight applied to tagging
data. In the case where the Negative Binomial is invoked, the following
expression is used:

$$\begin{matrix}
{L\left( {ObsRecap}_{r,y,a,s}^{k} \right) =} \\
{\lambda_{\text{Tagging}}\frac{\Gamma\left( {ObsRecap}_{r,y,a,s}^{k} + \eta \right)}{\Gamma(\eta)\Gamma\left( {ObsRecap}_{r,y,a,s}^{k} + 1 \right)}\left( \frac{\eta}{{Recap}_{r,y,a,s}^{k} + \eta} \right)^{\eta}\left( \frac{{Recap}_{r,y,a,s}^{k}}{{Recap}_{r,y,a,s}^{k} + \eta} \right)^{{ObsRecap}_{r,y,a,s}^{k}}} \\

\end{matrix}$$

Here, $\eta$ represents the estimated overdispersion parameter for
tagging data.

Under release conditioned dynamics, both recaptured and non-recaptured
states are fit to. Proportions of observed (${PObsRecap}_{r,y,a,s}^{k}$)
and expected recaptured (${PRecap}_{r,y,a,s}^{k})$individuals are given
by:

$$\begin{matrix}
{{PObsRecap}_{r,y,a,s}^{k} = \frac{{ObsRecap}_{r,y,a,s}^{k}}{{InitTag}^{k}}} \\

\end{matrix}$$

$${PRecap}_{r,y,a,s}^{k} = \frac{{Recap}_{r,y,a,s}^{k}}{{InitTag}^{k}}$$

${InitTag}^{k}$ denotes the total tags released for a given tag cohort
(combination of release region and year). Non-recaptured states can then
be written as:

$$\begin{matrix}
{{PObsNonRecap}_{r,y,a,s}^{k} = 1 - \sum\limits_{r}^{}{\sum\limits_{y}^{}{\sum\limits_{a}^{}{\sum\limits_{s}^{}{PObsRecap}_{r,y,a,s}^{k}}}}} \\

\end{matrix}$$

$${PNonRecap}_{r,y,a,s}^{k} = 1 - \sum\limits_{r}^{}{\sum\limits_{y}^{}{\sum\limits_{a}^{}{\sum\limits_{s}^{}{PRecap}_{r,y,a,s}^{k}}}}$$

where ${PObsNonRecap}_{r,y,a,s}^{k}$ and ${PNonRecap}_{r,y,a,s}^{k}$ are
the observed and expected non-recaptured states, respectively. These
states are then combined into a single vector of observed and expected
values:

$$\begin{matrix}
{\mathbf{O}_{Tagging}^{k} = \left\{ {\mathbf{P}\mathbf{O}\mathbf{b}\mathbf{s}\mathbf{R}\mathbf{e}\mathbf{c}\mathbf{a}}\mathbf{p}^{k},{PObsNonRecap}^{k} \right\}} \\

\end{matrix}$$

$$\mathbf{E}_{Tagging}^{k} = \{{\mathbf{P}\mathbf{R}\mathbf{e}\mathbf{c}\mathbf{a}}\mathbf{p}^{k},{PNonRecap}^{k}\}$$

If a Multinomial likelihood is assumed for release conditioned dynamics,
this is given by:

$$\begin{matrix}
{L\left( {ObsRecap}^{k} \right) = {\lambda_{\text{Tagging}}InitTag}^{k}\prod\limits_{i}^{}\left( E_{i,Tagging}^{k} \right)^{O_{i,Tagging}^{k}}} \\

\end{matrix}$$

Here, the subscript $i$ is used to generically denote a given element.
If a Dirichlet-multinomial with released-condition dynamics was assumed,
the tagging likelihood would be written as:

$$\begin{matrix}
{L\left( {ObsRecap}^{k} \right) =} \\
{\lambda_{\text{Tagging}} \cdot \frac{\Gamma\left( {InitTag}^{k}\  + 1 \right)}{{InitTag}^{k}\  \cdot O_{i,Tagging}^{k} + 1} \cdot \frac{\Gamma\left( {InitTag}^{k}\ \eta \right)}{\text{Γ(}{InitTag}^{k}\text{+}{InitTag}^{k}\eta\text{)}} \cdot} \\
{\prod\limits_{i}^{}\frac{\Gamma\left( {InitTag}^{k} \cdot O_{i,Tagging}^{k} + {InitTag}^{k} \cdot \eta \cdot E_{i,Tagging}^{k} \right)}{{InitTag}^{k} \cdot \eta \cdot E_{i,Tagging}^{k}}} \\

\end{matrix}$$

The $\eta$ parameter in the Dirichlet-multinomial likelihood represents
the overdispersion parameter for tagging data,

Under recapture-conditioned dynamics, tag shedding, tag induced
mortality, and tag reporting rates are assumed to be spatially-invariant
and do not need to be estimated, given that these terms cancel out in
the denominator (McGarvey and Feenstra, 2002). Unlike
release-conditioned dynamics, assuming recaptured-conditioned processes
does not require fitting to non-recaptured states. Thus, the observed
and expected recaptured proportions can be written as:

$$\begin{matrix}
{{PObsRecap}_{r,y,a,s}^{k} = \frac{{ObsRecap}_{r,y,a,s}^{k}}{\sum\limits_{r}^{}{\sum\limits_{a}^{}{\sum\limits_{s}^{}{ObsRecap}_{r,y,a,s}^{k}}}}} \\

\end{matrix}{{PRecap}_{r,y,a,s}^{k} = \frac{{Recap}_{r,y,a,s}^{k}}{\sum\limits_{r}^{}{\sum\limits_{a}^{}{\sum\limits_{s}^{}{Recap}_{r,y,a,s}^{k}}}}}{{TotRecap}_{y}^{k}\sum\limits_{r}^{}{\sum\limits_{a}^{}{\sum\limits_{s}^{}{ObsRecap}_{r,y,a,s}^{k}}}}$$

where recapture probabilities are normalized by the total number of
recaptures in a given year (${TotRecap}_{y}^{k}$).

### Parameter Priors and Process Error Penalties

#### Parameter Priors

Considering the complexity of integrated population models, several
priors can be specified to help inform the estimation of parameters by
providing additional knowledge. Beyond that, priors can also be
specified to aid in model stabilization (i.e., regularizing priors).
Priors can currently be specified for natural mortality, fishery and
survey catchability, fishery and survey selectivity, steepness,
recruitment proportions, movement rates, and tag reporting rates.

##### Natural Mortality

In the case of natural mortality, a lognormal prior is utilized:

$$\begin{matrix}
{P\left( \log\left( \text{NatMort}_{r,y,a,s} \right) \right) = \frac{1}{\sqrt{2\pi\sigma_{r,y,a,s}^{2\ {(Natmort)}}}\,}\exp\left( - \frac{\left\lbrack \log\left( \text{NatMort}_{r,y,a,s} \right) - log\left( \mu_{r,y,a,s}^{\text{(NatMort)}} \right) \right\rbrack^{2}}{2\sigma_{\text{NatMort}}^{2}\sigma_{r,y,a,s}^{2\ {(Natmort)}}} \right)} \\

\end{matrix}$$

where the variance of the prior is given by
$\sigma_{r,y,a,s}^{2\ {(Natmort)}}$, and
$\mu_{r,y,a,s}^{\text{(NatMort)}}$ denotes the prior mean.

##### Fishery and Survey Catchability

For fishery and survey catchability, a lognormal prior can also be
specified:

$$\begin{matrix}
{P\left( \log\text{(}q_{r,y,j}\text{)} \right) = \frac{1}{\sqrt{2\pi\sigma_{r,y,j}^{2\ {(qPrior)}}}\,}\exp\left( - \frac{\left\lbrack \log\left( q_{r,y,j} \right) - log\left( \mu_{r,y,j}^{\text{(qPrior)}} \right) \right\rbrack^{2}}{2\sigma_{\text{NatMort}}^{2}\sigma_{r,y,j}^{2\ {(qPrior)}}} \right)} \\

\end{matrix}$$

where $j$ here represents either a fishery or survey fleet,
$\sigma_{r,y,j}^{2\ {(qPrior)}}$ is the variance of the prior, and
$\mu_{r,y,j}^{\text{(qPrior)}}$ indicates the prior mean for
catchability.

##### Fishery and Survey Selectivity

In general, selectivity priors can be utilized to serve as regularizing
priors to facilitate stable parameter estimation (Monnahan, 2024). These
priors are assumed to be lognormal and are applied to the selectivity
parameters themselves:

$$\begin{matrix}
{P\left( \log\left( \theta_{r,p,y,s,j} \right) \right) = \frac{1}{\sqrt{2\pi{\sigma^{2}}_{r,p,y,s,j}^{\text{(Sel)}}}\,}\exp\left( - \frac{\left\lbrack \ln\left( \theta_{r,p,y,s,j} \right) - ln\left( \mu_{r,p,y,s,j}^{\text{(Sel)}} \right) \right\rbrack^{2}}{2{\sigma^{2}}_{r,p,y,s,j}^{\text{(Sel)}}} \right)} \\

\end{matrix}$$

where $\theta_{r,p,y,s,j}$ is a selectivity parameter for a given
functional form specified, ${\sigma^{2}}_{r,p,y,s,j}^{\text{(Sel)}}$ is
the prior variance, and $\mu_{r,p,y,s,j}^{\text{(Sel)}}$ is the prior
mean for the specific selectivity parameter.

##### Steepness

If a Beverton-Holt stock recruitment relationship is assumed, priors for
steepness can be specified. Currently, a scaled beta prior (bounded
between 0.2 and 1) can be invoked:

$$\begin{matrix}
{a_{r}^{(h)} = \left( \frac{\mu_{r}^{\text{(}\text{h}\text{)}} - 0.2}{1 - 0.2} \right)\left( \frac{\sigma_{r}^{(h)}}{1 - 0.2} \right)^{2}} \\

\end{matrix}{b_{r}^{(h)} = \left\lbrack 1 - \ \left( \frac{\mu_{r}^{\text{(}\text{h}\text{)}} - 0.2}{1 - 0.2} \right) \right\rbrack\left\lbrack \left( \frac{\sigma_{r}^{(h)}}{1 - 0.2} \right)^{2} \right\rbrack}{P\left( h_{r} \right) = \frac{\Gamma\left( a_{r}^{(h)} \right)\Gamma\left( b_{r}^{(h)} \right)}{\Gamma\left( {a_{r}^{(h)}}_{r} + b_{r}^{(h)} \right)}h_{r}^{a - 1}\left( 1 - h_{r} \right)^{b - 1}}$$

Here, $a_{r}^{(h)}$ and $b_{r}^{(h)}$ are parameters of the beta
distribution, $\mu_{r}^{\text{(}\text{h}\text{)}}$ is the prior mean
steepness for a given region (bounded between 0.2 and 1) while
$\sigma_{r}^{(h)}$ is the standard deviation for these priors.

##### Recruitment Proportions

Regional recruitment is derived by apportioning a global recruitment
parameter using regional recruitment proportions (i.e.,
$\mu^{\text{Rec}} \cdot \zeta_{r}$). Here, $\zeta_{r}$ is derived via a
multinomial logit transformation and Dirichlet priors can be used to
help constrain estimation:

$$\begin{matrix}
{P({\mathbf{ζ}}) = \frac{\Gamma\left( \sum\limits_{r}^{n_{r}}\alpha_{r} \right)}{\prod\limits_{r}^{n_{r}}{\Gamma\left( \alpha_{r} \right)}}\prod\limits_{r = 1}^{n_{r}}\zeta_{r}^{\alpha_{r} - 1}} \\

\end{matrix}$$

${\mathbf{ζ}} = \{\zeta_{1},\zeta_{2},\ldots,\zeta_{n_{r}}\}$ are the
estimated recruitment proportions across regions, and $\alpha_{r}$ is
the concentration parameter governing the spread of the Dirichlet
distribution.

##### Movement

Likewise, priors on movement values can be assumed to arise from a
Dirichlet process:

$$\begin{matrix}
{P\left( \mathbf{M}_{.,k,y,a,s} \right) = \frac{\Gamma\left( \sum\limits_{r}^{n_{r}}c_{r,k,y,a,s} \right)}{\prod\limits_{r}^{n_{r}}{\Gamma\left( c_{r,k,y,a,s} \right)}}\prod\limits_{r = 1}^{n_{r}}M_{r,k,y,a,s}^{c_{r,k,y,a,s} - 1}} \\

\end{matrix}$$

where
$\mathbf{M}_{.,k,y,a,s}\ \epsilon\ \{ M_{1,k,y,a,s,\ \ldots,\ }M_{n_{r,k,y,a,s}}\}$,
$r$ is the origin region, $k$ is the destination, and $c_{r,k,y,a,s}$
are the concentration parameters that control the Dirichlet
distribution.

##### Tag Reporting Rates

Two types of priors can be specified for tag reporting rates. In
particular, a symmetric beta distribution is applied:

$$\begin{matrix}
{P\left( \beta_{r,y} \right) = \left( \beta_{r,y} + 1e - 4 \right)^{\sigma_{r,y}^{(\beta)}}\left( 1 - \beta_{r,y} + 1e - 4 \right)^{\sigma_{r,y}^{(\beta)}}} \\

\end{matrix}$$

Here, $\sigma_{r,y}^{(\beta)}$ determines the scale of the tag reporting
parameter and determines how strongly to penalize estimates when they
approach the bounds of $\lbrack 0,1\rbrack$. Smaller values of
$\sigma_{r,y}^{(\beta)}$ result in larger penalties, and vice versa.
Thus, estimates of tag reporting rates are penalized if they are either
close to 0 or 1.

Tag reporting rate priors can also be specified as a standard beta
distribution:

$$\begin{matrix}
{a_{r}^{(\beta)} = \frac{\mu_{r,y}^{(\beta)}}{\left( \sigma_{r,y}^{(\beta)} \right)^{2}}} \\

\end{matrix}{b_{r}^{(\beta)} = \left\lbrack \frac{1 - \mu_{r,y}^{(\beta)}}{\left( \sigma_{r,y}^{(\beta)} \right)^{2}} \right\rbrack}{P\left( \beta_{r,y} \right) = \frac{\Gamma\left( a_{r}^{(\beta)} \right)\Gamma\left( b_{r}^{(\beta)} \right)}{\Gamma\left( a_{r}^{(\beta)} + b_{r}^{(\beta)} \right)}\beta_{r,y}^{a - 1}\left( 1 - \beta_{r,y} \right)^{b - 1}}$$

Here, $a_{r}^{(h)}$ and $b_{r}^{(h)}$ are the parameters describing the
beta distribution, $\mu_{r,y}^{(\beta)}$ denotes the prior mean for tag
reporting rates and $\sigma_{r,y}^{(\beta)}$ is the prior standard
deviation.

#### Process Error Penalties

In addition to priors, penalties are also utilized to aid in the
estimation of process errors (either penalized likelihood or integrating
random effects via Laplace Approximation are possible). Currently,
process errors can be specified to arise for initial age deviations,
recruitment, fishing mortality, fishery and survey selectivity, and
movement.

##### Initial Age Deviations

To estimate non-equilibrium initial age deviations, multiplicative
deviations can be specified:

$$L\left( \epsilon_{r,i}^{\text{Init}} \right) = \frac{1}{\sqrt{2\pi\sigma_{\text{Init}}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,i}^{\text{Init}} \right)^{2}}{2\sigma_{\text{Init}}^{2}} \right)$$

where deviations arise from a normal distribution with a mean of 0 and
variance of $\sigma_{\text{Init}}^{2}$.

##### Recruitment Deviations

Annual recruitment deviations can also be specified, where
multiplicative deviations are assumed:

$$L\left( \epsilon_{r,y}^{\text{Rec}} \right) = \frac{1}{\sqrt{2\pi}\sigma_{\text{Rec}}^{2}\,}\exp\left( - \frac{\left( \epsilon_{r,y}^{\text{Rec}} \right)^{2}}{2\sigma_{\text{Rec}}^{2}} \right)$$

and deviations are assumed to normally distributed, with a mean of 0 and
variance of $\sigma_{\text{Rec}}^{2}$.

##### Fishing Mortality Deviations

Fishing mortality deviations similarly assumes multiplicative
deviations:

$$L\left( \epsilon_{r,y,f}^{\text{Fsh}} \right) = \frac{1}{\sqrt{2\pi{\sigma_{r,f}^{2}}_{\text{Fsh}}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,f}^{\text{Fsh}} \right)^{2}}{2{\sigma_{r,f}^{2}}_{\text{Fsh}}} \right)$$

Here, fishing mortality deviations are assumed to arise from a normal
distribution, with mean 0 and a variance of
${\sigma_{r,f}^{2}}_{\text{Fsh}}$.

##### Fishery and Survey Selectivity

A variety of process error parameterizations can be specified for
fishery and survey selectivity. Across all parameterizations,
multiplicative deviations are assumed. In the most basic case, iid
deviations can be assumed to vary about a parameter on a given
selectivity functional form:

$$L\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} \right) = \frac{1}{\sqrt{2\pi\sigma_{r,i,s,j,Sel}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} \right)^{2}}{2\sigma_{r,i,s,j,Sel}^{2}} \right)$$

where $\epsilon_{r,y,i,s,j}^{\text{Sel}}$ are selectivity deviations
about a given parameter for region $r$, year $y$, parameter $i$, sex
$s$, and fleet $j$. Deviations are assumed to have a mean of 0 and a
variance of $\sigma_{r,i,s,j,Sel}^{2}$, constrained by a normal
distribution.

Extending the iid case, random walk selectivity deviations can also be
specified about a given parameter, assuming a normal distribution:

$$L\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} \right) = \left\{ \begin{matrix}
{\frac{1}{\sqrt{2\pi \cdot 5^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y = 1,i,s,j}^{\text{Sel}} \right)^{2}}{2{\cdot 5}^{2}} \right),\ ify = 1)} \\
{\frac{1}{\sqrt{2\pi\sigma_{r,i,s,j,Sel}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} - \epsilon_{r,y - 1,i,s,j}^{\text{Sel}} \right)^{2}}{2\sigma_{r,i,s,j,Sel}^{2}} \right),ify > 1)} \\

\end{matrix} \right.\ $$

where process error deviations for the first year are initialized with a
large variance. Following the first year, process error deviations
follow a random walk process with a mean conditional on the previous
year’s value ($\epsilon_{r,y - 1,i,s,j}^{\text{Sel}}$) and a variance of
$\sigma_{r,i,s,j,Sel}^{2}$.

In addition to being constrained by a normal distribution, both iid and
random walk cases have an optional additional smoothness penalty
applied:

$$P_{SelYrSmooth} = \sum\limits_{r = 1}^{n_{r}}{\sum\limits_{j = 1}^{n_{j}}{\sum\limits_{s = 1}^{n_{s}}{\sum\limits_{y = 2}^{n_{y}}{\sum\limits_{i = 1}^{n_{i}}\left( \epsilon_{r,y,i,s,j}^{\text{Sel}} - \epsilon_{r,y - 1,i,s,j}^{\text{Sel}} \right)^{2}}}}}$$

which serves to ensure selectivity values do not vary drastically from
year-to-year.

Additionally, semi-parametric deviations can also be specified. In
total, there are three options that can be utilized, two of which allow
age, year, and cohort correlations, while one allows for only age or
length and year correlations. In the case where age, year, and cohort
correlations are specified (note that this is only possible when
age-based selectivity is specified), marginal stationary variance and a
conditional non-stationary variance can be invoked. The primary
difference between the two is that the marginal variance version does
not have a closed form solution for $\mathbf{\Omega}$ (see equations
below) and needs to be iteratively solved, while the conditional
variance version has a closed form solution. The following equations
describe the conditional variance version:

$${\mathbf{ϵ}}_{r,s,j}^{\text{Sel}} = vec\left( \epsilon_{r,y,a,s,j}^{\text{Sel}} \right)$$

where we vectorize the selectivity deviations across its year and age
dimensions. These deviations are then assumed to arise from a
multivariate normal distribution (or Gaussian Markov Random Field) with
a covariance matrix ($\mathbf{\Sigma} = \mathbf{Q}^{- 1}$) determined
by:

$$\begin{matrix}
{\text{Q} = \left( \text{I} - \left( \text{B} \right)^{T} \right)\mathbf{\Omega}\left( \text{I} - \text{B} \right)} \\
{\text{diag}(\mathbf{\Omega}) = \sigma_{r,s,j,Sel}^{- 2}} \\

\end{matrix}$$

Here, $\text{I}$ is an identity matrix and $\mathbf{\Omega}$ is a
diagonal matrix that determines the variance of the multivariate normal
process. $\text{B}$ is a square matrix representing the partial effect
of ${\mathbf{ϵ}}_{r,s,j}^{\text{Sel}}$ on preceding ages and/or years,
governed by partial correlation coefficients for ages, years, and
cohorts. To demonstrate the formulation of $\text{B}$, a simplified
example is provided with rows representing ages $a\ \epsilon\ \{ 1,2\}$
and columns representing years $y\ \epsilon\ \{ 1,2\}$. In this example,
$\text{B}$ is a $4\ x\ 4$ matrix, where both the rows and columns
represent combinations of age and year. For instance, $B_{1,2}$ captures
the correlation within the same year between age 1 in year 1, and age 2
in year 1, whereas $B_{4,1}$ constructs the correlation within the same
cohort between age 1 in year 1, and age 2 in year 2:

$$\begin{matrix}
{\mathbf{B} = \begin{bmatrix}
1 & 0 & 0 & 0 \\
\rho_{y} & 0 & 0 & 0 \\
\rho_{a} & 0 & 0 & 0 \\
\rho_{c} & \rho_{a} & \rho_{y} & 0 \\
 & & & 
\end{bmatrix}} \\

\end{matrix}$$

where $\rho_{y}$, $\rho_{a}$, and $\rho_{c}$ are parameters describing
the partial autocorrelation among years within a given age, among ages
within a given year, and years within a cohort, respectively. The
multivariate likelihood is then defined as:

$$L\left( {\mathbf{ϵ}}_{r,s,j}^{\text{Sel}} \right) = \frac{|\mathbf{Q}|^{1/2}}{(2\pi)^{\frac{n}{2}}}\exp\left( - \frac{1}{2}\left( {\mathbf{ϵ}}_{r,s,j}^{\text{Sel}} \right)^{T}\mathbf{Q}\left( {\mathbf{ϵ}}_{r,s,j}^{\text{Sel}} \right) \right)$$

If age or length and year correlations are specified (i.e., a
two-dimensional autoregressive structure), a multivariate normal
likelihood is similarly assumed (i.e., the equation described above),
but the covariance structured of this process is defined as:

$$\text{Q}^{- 1} = \frac{\sigma_{r,s,j,Sel}^{2}}{\left( 1 - \rho_{y} \right)^{2}\left( 1 - \rho_{b} \right)^{2}}\text{R}_{y} \otimes \text{R}_{b}$$

where $\rho_{y}$ and $\rho_{b}$ are correlation coefficients across
years and bins, respectively. $\text{R}_{y}$ is a $n_{y} \times n_{y}$
matrix representing a first-order autoregressive process across years
and $\text{R}_{b}$ is a $n_{b} \times n_{b}$ matrix representing a
first-order autoregressive process across bins.

Moreover, when semi-parametric deviations are specified, additional
optional penalties can be applied across bins and years to enforce
curvature control:

$$P_{SelBinCurve} = \sum\limits_{r = 1}^{n_{r}}{\sum\limits_{j = 1}^{n_{j}}{\sum\limits_{s = 1}^{n_{s}}{\sum\limits_{y = 1}^{n_{y}}{\sum\limits_{b = 2}^{n_{b} - 1}\left( \log\left( {Sel}_{r,y,b + 1,s,j} \right) - 2\log\left( {Sel}_{r,y,b,s,j} \right)\log\left( {Sel}_{r,y,b - 1,s,j} \right) \right)^{2}}}}}$$

$$P_{SelYrCurve} = \sum\limits_{r = 1}^{n_{r}}{\sum\limits_{j = 1}^{n_{j}}{\sum\limits_{s = 1}^{n_{s}}{\sum\limits_{y = 2}^{n_{y} - 1}{\sum\limits_{b = 1}^{n_{b}}\left( \log\left( {Sel}_{r,y + 1,b,s,j} \right) - 2\log\left( {Sel}_{r,y,b,s,j} \right)\log\left( {Sel}_{r,y - 1,b,s,j} \right) \right)^{2}}}}}$$

##### Movement

Process error deviations can also be specified for movement parameters.
Currently, only iid deviations are allowed, where deviations are placed
upon fixed-effects movement parameters in logit space:

$$\omega_{r,k,y,a,s} = \omega_{r,k,s} + \epsilon_{r,k,y,a,s}^{\text{Move}}$$

where $\omega_{r,k,y,a,s}$ are movement parameters in logit space for
origin $r$, destination $k$, year $y$, age $a$, and sex $s$, and
$\omega_{r,k,s}$ are the corresponding fixed effects (the mean logits).
Here, $\epsilon_{r,k,y,a,s}^{\text{Move}}$ are movement deviations
specified in logit space. Deviations then arise from a normal
distribution with process variation constrained by
$\sigma_{r,a,s,Move}^{2}$:

$$L\left( \epsilon_{r, - 1,y,a,s}^{\text{Move}} \right) = \frac{1}{\sqrt{2\pi\sigma_{r,a,s,Move}^{2}}\,}\exp\left( - \frac{\left( \epsilon_{r, - 1,y,a,s}^{\text{Move}} \right)^{2}}{2\sigma_{r,a,s,Move}^{2}} \right)$$

Thus, process variation can be distinct for a given origin region, age,
and sex. Because logits are parameterized with a reference category,
only $n_{r} - 1$ logits (and their deviations) are estimated per origin
$r$, ensuring the sum-to-one constraint on destination probabilities.

#### Joint Likelihood

Lastly, the joint likelihood to be minimized represents the sum of all
observational likelihood components, priors, and penalties defined
above:

$$\begin{matrix}
{\text{Joint Likelihood} = \sum\text{Observation Likelihoods} + \sum\text{Priors} + \sum\text{Penalties}} \\

\end{matrix}$$

Note that some of these components may be zero (i.e., if no priors are
used) depending on the configuration of the model.

## References

Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H., Bell, B., 2016. TMB:
Automatic Differentiation and Laplace Approximation. J. Stat. Soft. 70.
<https://doi.org/10.18637/jss.v070.i05>

Mace, P.M., Doonan, I.J., 1988. A Generalised Bioeconomic Simulation
Model for Fish Population Dynamics. MAFFish, N.Z. Ministry of
Agriculture and Fisheries.

McGarvey, R., Feenstra, J.E., 2002. Estimating rates of fish movement
from tag recoveries: conditioning by recapture. Can. J. Fish. Aquat.
Sci. 59, 1054–1064. <https://doi.org/10.1139/f02-080>

Methot, R.D., Taylor, I.G., 2011. Adjusting for bias due to variability
of estimated recruitments in fishery assessment models. Can. J. Fish.
Aquat. Sci. 68, 1744–1760. <https://doi.org/10.1139/f2011-092>

Monnahan, C.C., 2024. Toward good practices for Bayesian data-rich
fisheries stock assessments using a modern statistical workflow.
Fisheries Research 275, 107024.
<https://doi.org/10.1016/j.fishres.2024.107024>

Thorson, J.T., Johnson, K.F., Methot, R.D., Taylor, I.G., 2017.
Model-based estimates of effective sample size in stock assessment
models using the Dirichlet-multinomial distribution. Fisheries Research
192, 84–93. <https://doi.org/10.1016/j.fishres.2016.06.005>
