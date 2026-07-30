# Initialize Numbers-at-Age (NAA) for a Population Model

Computes equilibrium initial numbers-at-age (NAA) for a structured
population model across populations, regions, sexes, and ages. Several
initialization methods are supported, ranging from numerical iteration
to analytical geometric-series solutions that optionally incorporate
seasonal movement.

## Usage

``` r
Get_Init_NAA(
  init_age_strc,
  init_iter,
  n_regions,
  n_pop,
  n_sexes,
  n_ages,
  n_seas,
  n_fish_fleets,
  seasdur,
  rec_seas_prop,
  natmort,
  init_F,
  dmr,
  fish_sel,
  ret_sel,
  R0_r,
  sexratio,
  Movement,
  do_recruits_move,
  ln_InitDevs,
  Mrate = NULL,
  move_timing = 0
)
```

## Arguments

- init_age_strc:

  Integer specifying the initialization method:

  - `0` Iterative equilibrium solution

  - `1` Scalar geometric-series solution (no movement in any age)

  - `2` Matrix geometric-series solution (movement allowed)

  - `3` Hybrid solution: movement in ages \< plus group, scalar solution
    for plus group

- init_iter:

  Integer; number of annual iterations used when `init_age_strc = 0`.

- n_regions:

  Integer; number of spatial regions.

- n_pop:

  Integer; number of populations.

- n_sexes:

  Integer; number of sexes.

- n_ages:

  Integer; number of age classes (including the plus group).

- n_seas:

  Integer; number of seasons per year.

- n_fish_fleets:

  Integer; number of fishing fleets.

- seasdur:

  Numeric vector (`n_seas`) giving the fraction of the year represented
  by each season.

- rec_seas_prop:

  Matrix (`n_pop x n_seas`) giving seasonal recruitment proportions.

- natmort:

  Array (`n_pop x n_regions x n_ages x n_sexes`) of natural mortality
  rates.

- init_F:

  Numeric array (`n_regions x n_seas x n_fish_fleets`) giving fishing
  mortality applied in each region, season, and fleet during
  initialization. Set to zero for an unfished population.

- dmr:

  Numeric array (`n_regions x n_seas x n_fish_fleets`) giving discard
  mortality rate applied in each region, season, and fleet during
  initialization (first year). Set to zero for an unfished population.

- fish_sel:

  Array
  (`n_pop x n_regions x n_seas x n_ages x n_sexes x n_fish_fleets`) of
  total fishery selectivity at age.

- ret_sel:

  Array
  (`n_pop x n_regions x n_seas x n_ages x n_sexes x n_fish_fleets`) of
  retained fishery selectivity at age.

- R0_r:

  Matrix (`n_pop x n_regions`) giving unfished recruitment allocated to
  each region.

- sexratio:

  Array (`n_pop x n_regions x n_sexes`) giving the proportion of
  recruits by sex.

- Movement:

  Array (`n_pop x origin x destination x n_seas x n_ages x n_sexes`)
  containing seasonal movement probabilities.

- do_recruits_move:

  Integer indicator:

  - `0` Recruits do not move during their first year

  - `1` Recruits move according to the movement matrix

- ln_InitDevs:

  Array (`n_pop x n_regions x (n_ages - 1)`) containing log-scale
  deviations applied to ages 2 through \\A\\.

## Details

The resulting age structure represents an equilibrium population under
constant recruitment (\\R_0\\), mortality, fishing mortality, and
movement.

Initial numbers-at-age are derived assuming constant recruitment
(\\R_0\\) and constant mortality and movement.

Let

- \\N\_{p,r,a,s}\\ denote numbers-at-age

- \\M\_{p,r,a}\\ denote natural mortality

- \\F\_{r,a,s}\\ denote total fishing mortality (retained + dead
  discards)

- \\Z = M + F\\ denote total mortality

Recruitment at age 1 is

\$\$ N\_{p,r,1,s} = R\_{0,p,r} \times sexratio\_{p,r,s} \$\$

Within-season survival follows

\$\$ N\_{p,r,a,s+1} = N\_{p,r,a,s}\exp(-Z\_{p,r,a,s}) \$\$

Ages advance at the end of the final season of the year:

\$\$ N\_{p,r,a+1,1} = N\_{p,r,a,n\_{seas}} \exp(-Z\_{p,r,a,n\_{seas}})
\$\$

The plus group accumulates survivors from the terminal age:

\$\$ N\_{A^+} = N\_{A-1} e^{-Z\_{A-1}} + N\_{A^+} e^{-Z\_{A^+}} \$\$

Fishing mortality at age is decomposed into retained and dead discard
components, summed across all fleets:

\$\$ F\_{p,r,s,a} = \sum\_{f=1}^{n_f} F^{init}\_{r,s,f} \left\[
sel\_{p,r,s,a,f} \cdot ret\_{p,r,s,a,f} + sel\_{p,r,s,a,f} \cdot (1 -
ret\_{p,r,s,a,f}) \cdot dmr\_{r,s,f} \right\] \$\$

where \\sel\\ is total fishery selectivity, \\ret\\ is retention
selectivity, and \\dmr\\ is the discard mortality rate.

\### Scalar geometric-series solution

When movement is absent, equilibrium abundance follows

\$\$ N_a = N_1 \exp\left(-\sum\_{i=1}^{a-1} Z_i \right) \$\$

The plus group has a closed-form solution

\$\$ N\_{A^+} = \frac{N\_{A-1} e^{-Z\_{A-1}}} {1 - e^{-Z\_{A^+}}} \$\$

\### Matrix geometric-series solution

When movement occurs, survival and movement are combined into seasonal
transition matrices:

\$\$ \mathbf{T}\_a = \prod\_{s=1}^{n\_{seas}}
\mathbf{M}\_{a,s}\mathbf{S}\_{a,s} \$\$

where

- \\\mathbf{M}\\ is the movement transition matrix

- \\\mathbf{S}\\ is a diagonal matrix of survival probabilities

The plus group equilibrium satisfies

\$\$ \mathbf{N}\_{A^+} = (\mathbf{I} - \mathbf{T}\_{A^+})^{-1}
\mathbf{T}\_{A-1} \mathbf{N}\_{A-1} \$\$

The iterative method (`init_age_strc = 0`) numerically applies the full
seasonal population dynamics repeatedly until the population converges
to equilibrium.

After equilibrium is derived, log-scale initial age deviations
(`ln_InitDevs`) are applied multiplicatively to ages \\2,\dots,A\\.

Initial numbers-at-age are derived assuming constant recruitment
(\\R_0\\) and constant mortality and movement.

Let

- \\N\_{p,r,a,s}\\ denote numbers-at-age

- \\M\_{p,r,a}\\ denote natural mortality

- \\F\_{r,a,s}\\ denote total fishing mortality (retained + dead
  discards)

- \\Z = M + F\\ denote total mortality

Recruitment at age 1 is

\$\$ N\_{p,r,1,s} = R\_{0,p,r} \times sexratio\_{p,r,s} \$\$

Within-season survival follows

\$\$ N\_{p,r,a,s+1} = N\_{p,r,a,s}\exp(-Z\_{p,r,a,s}) \$\$

Ages advance at the end of the final season of the year:

\$\$ N\_{p,r,a+1,1} = N\_{p,r,a,n\_{seas}} \exp(-Z\_{p,r,a,n\_{seas}})
\$\$

The plus group accumulates survivors from the terminal age:

\$\$ N\_{A^+} = N\_{A-1} e^{-Z\_{A-1}} + N\_{A^+} e^{-Z\_{A^+}} \$\$

Fishing mortality at age is decomposed into retained and dead discard
components, summed across all fleets:

\$\$ F\_{p,r,s,a} = \sum\_{f=1}^{n_f} F^{init}\_{r,s,f} \left\[
sel\_{p,r,s,a,f} \cdot ret\_{p,r,s,a,f} + sel\_{p,r,s,a,f} \cdot (1 -
ret\_{p,r,s,a,f}) \cdot dmr\_{r,s,f} \right\] \$\$

where \\sel\\ is total fishery selectivity, \\ret\\ is retention
selectivity, and \\dmr\\ is the discard mortality rate.

\### Scalar geometric-series solution

When movement is absent, equilibrium abundance follows

\$\$ N_a = N_1 \exp\left(-\sum\_{i=1}^{a-1} Z_i \right) \$\$

The plus group has a closed-form solution

\$\$ N\_{A^+} = \frac{N\_{A-1} e^{-Z\_{A-1}}} {1 - e^{-Z\_{A^+}}} \$\$

\### Matrix geometric-series solution

When movement occurs, survival and movement are combined into seasonal
transition matrices:

\$\$ \mathbf{T}\_a = \prod\_{s=1}^{n\_{seas}}
\mathbf{M}\_{a,s}\mathbf{S}\_{a,s} \$\$

where

- \\\mathbf{M}\\ is the movement transition matrix

- \\\mathbf{S}\\ is a diagonal matrix of survival probabilities

The plus group equilibrium satisfies

\$\$ \mathbf{N}\_{A^+} = (\mathbf{I} - \mathbf{T}\_{A^+})^{-1}
\mathbf{T}\_{A-1} \mathbf{N}\_{A-1} \$\$

The iterative method (`init_age_strc = 0`) numerically applies the full
seasonal population dynamics repeatedly until the population converges
to equilibrium.

After equilibrium is derived, log-scale initial age deviations
(`ln_InitDevs`) are applied multiplicatively to ages \\2,\dots,A\\.
