# EBS Pacific cod bridge data (Model 24.1)

The 2024 eastern Bering Sea Pacific cod assessment (Model 24.1: one area
and one sex, Richards growth with annual deviations on the length at the
young reference age and on K carried cohort by cohort, length-based
double normal selectivity with two fishery blocks and annual deviations
on the survey width, compositions recorded on coarser length bins than
the population carries, two ageing error definitions) parsed from its
input files, with the assessment's report quantities attached under
`$ss3`, its year-by-year weight and fecundity at age under `$wtatage`,
and its maximum likelihood estimates under `$mle`. Built by
`dev/pcod_bridge/R/build_pcod_data.R` from the AFSC EBS_PCOD repository.
Used by `tests/testthat/helper-bridge_ebs_pcod.R` and the Pacific cod
regression test, which evaluates SPoRC at the assessment's own estimate.

## Usage

``` r
sgl_rg_ebs_pcod_data
```

## Format

A named list of observation arrays, biological inputs, the assessment's
configuration, `$LenBinMap` (the population-to-data length bin map),
`$MatAA` (maturity at age as the share of the weight at age that
spawns), `$wtatage` (weight and fecundity at age by year), `$mle`
(parameter values) and `$ss3` (report quantities: numbers at age,
spawning biomass, recruitment, total biomass, catch, the index, growth
parameters by year, selectivity at length and age, expected
compositions, likelihood components and deviation penalties).

## Source

<https://github.com/afsc-assessments/EBS_PCOD>
