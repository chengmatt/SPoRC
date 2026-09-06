# Development Roadmap

This document records what `SPoRC` development is working towards, what
is saved for later to be revisited, and what is likely not in the
development plan (out of scope) .

Four tracks are likely to be developed in the future: numbers at age as
a state space process, close-kin mark-recapture, electronic tagging, and
process error on the remaining parameter blocks.

Work to date has involved the following released versions:

## Released versions

| Version | Date | What it added |
|----|----|----|
| 1.0.0 | 2025-11-24 | First public release. Age, sex, season and region structured dynamics, closed loop simulation and management strategy evaluation. |
| 1.1.0 | 2026-03-31 | Movement estimated as a continuous time Markov chain with preference functions. |
| 1.2.0.9000 | in development | Population specific (natal homing) and seasonal dynamics, discarding and retention, time varying and semi-parametric growth, conditional age at length, at age data sources, estimated index observation error, per fleet ageing error. |

## Numbers at age as a state space process

Recruitment is the only stochastic element of the population state. From
age 2 onwards, numbers at age are a deterministic consequence of
recruitment, mortality and movement, so any misspecification in natural
mortality, in selectivity or in the reported catch is absorbed by the
recruitment deviations. Estimating process error on the whole numbers at
age surface, as in SAM and WHAM, lets that misspecification enter at the
age it occurs instead. (This has now been implemented).

## Close-kin mark-recapture

Close-kin mark-recapture estimates absolute abundance from the rate at
which pairs of sampled fish turn out to be related. The quantities it
needs are already in the report file: numbers at age by population,
region, year and sex, maturity at age, and weight at age, from which
relative reproductive output follows. What is missing is the pairwise
comparison data and a likelihood over it.

## Electronic tagging

Conventional tags inform movement through recaptures. Recaptures are
sparse, they happen only where the fishery is, and they arrive filtered
through a reporting rate. Archival and satellite tags record position
through the fish’s life independently of the fishery, so they inform the
movement rates themselves rather than the movement rates confounded with
fishing processes.

## Process error on the remaining parameter blocks

Selectivity and growth already take a deviation surface with a choice of
covariance, and fishing mortality takes independent, random walk or
first order autoregressive deviations. Natural mortality and
catchability are constant within a time block, and regional recruitment
deviations are independent. Thus, another priority of `SPoRC` will be to
add additional process error functionality to these data sources.

## To be revisited

Environmental covariates on recruitment, natural mortality, growth and
catchability. Movement already accepts covariates through
`cont_vary_movement`, and extending the same interface to the other
processes is a natural addition. It is not a near term priority and will
be revisited.

## Outside current scope

| Capability | Rationale |
|----|----|
| Multispecies models and predation mortality | `SPoRC` is designed for single-species population dynamics, potentially with multiple populations or regions. Multispecies interactions and predation mortality are therefore outside the current scope of the framework. |
| Fully length-structured dynamics | `SPoRC` uses age as the primary population-structure dimension. Length-based observations and selectivity can be represented through age–length relationships without requiring the population dynamics themselves to be fully length structured. A fully length-structured implementation would substantially expand the model structure and computational burden without being necessary for the primary applications of `SPoRC`. |
| Hermaphroditism and sex change | `SPoRC` assumes sex is fixed after recruitment. Given that hermaphroditism and sex change are relevant to relatively few species and are not a primary focus of the intended applications, they are outside the current scope of the framework. |
