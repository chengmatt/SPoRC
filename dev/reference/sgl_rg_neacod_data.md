# Northeast Arctic cod data for the SAM case study

Inputs for the Northeast Arctic cod assessment fit with the State-space
Assessment Model (SAM), 1946 to 2024 at ages 3 to 15 with a plus group:
catch at age from one fishery and survey indices at age from five
surveys, with the estimates SAM returns for the same data, used to check
the bridge.

## Usage

``` r
sgl_rg_neacod_data
```

## Format

A list with components `inputs` (dimensions, weight, maturity and
natural mortality at age, the catch and survey observations at age with
their use flags, survey timing, and SAM's observation error,
catchability and F key matrices shifted to start at one), `sam` (SAM's
estimates of log numbers and log F at age, log catchability, the
observation and process standard deviations, the F correlation across
ages, spawning biomass and the negative log likelihood, under compound
symmetry in F) and `fit_summary` (SAM's negative log likelihood under
independent, compound symmetric and AR1 F increments).

## Source

ICES Arctic Fisheries Working Group input files for Northeast Arctic
cod, catches 1946 to 2023 and surveys to 2024, fit with `samjr`, the R
implementation of SAM (Nielsen, A., Berg, C. W. 2014. Estimation of
time-varying selectivity in stock assessments using state-space models.
Fisheries Research 158: 96-101).
