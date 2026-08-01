# Normal prior on global R0

Called once from the "Recruitment R0 (Prior)" section of `SPoRC_rtmb.R`.
Returns a scalar contribution added directly (unweighted by `Wt_Rec`)
into the joint negative log-likelihood, alongside the other scalar
priors (`M_nLL`, `h_nLL`, etc.); it is not part of the `Rec_nLL`
recruitment-deviation array, since it penalizes a single global
parameter rather than a per-year deviation.

## Usage

``` r
get_r0_prior(r0_prior, ln_global_R0)
```

## Arguments

- r0_prior:

  Data frame with columns `pop`, `mu` (prior mean R0, natural scale),
  `sd` (prior SD, log scale), one row per penalized population.

- ln_global_R0:

  Vector `[pop]` of log mean recruitment (R0).

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `r0_prior`.
