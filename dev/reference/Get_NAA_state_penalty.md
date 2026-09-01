# State-space numbers at age

Scores the realized innovation of the centered numbers-at-age state,
\\\eta = \log N - \log \hat{N}\\, where \\\hat{N}\\ is the deterministic
mortality and ageing prediction the dynamics computed into `NAA_pred`
before the state overwrote `NAA`.

## Usage

``` r
Get_NAA_state_penalty(
  ln_NAA,
  NAA_pred,
  sigmaNAA,
  naa_re_ages,
  naa_re_yrs,
  NAA_re = 1,
  NAA_pe_pars = NULL,
  NAA_re_region = 0,
  NAA_region_corr_pars = NULL,
  NAA_re_pop = 0,
  NAA_pop_corr_pars = NULL,
  NAA_re_sex = 0,
  NAA_sex_corr_pars = NULL
)
```

## Arguments

- ln_NAA:

  Array `[pop, region, year, age, sex]` of log numbers at age.

- NAA_pred:

  Array of the same shape holding the deterministic prediction.

- sigmaNAA:

  Array of the same shape holding the process error standard deviation
  for each cell, already expanded from its blocking structure.

- naa_re_ages:

  Integer vector of age indices the state is active over.

- naa_re_yrs:

  Integer vector of year indices the state is active over.

- NAA_re:

  Integer code for the structure over the age-year grid. `1`
  independent, `2` AR(1) over ages, `3` AR(1) over years with ages
  independent, `4` separable AR(1) over ages and years, `5` and `6` the
  three-dimensional Gaussian Markov random field on the conditional and
  the marginal variance respectively.

- NAA_pe_pars:

  Array `[pop, region, 3, sex]` of correlation parameters on the
  unconstrained scale, read as age, year and cohort. Unused under
  `NAA_re = 1`.

- NAA_re_region:

  Integer code for the structure across regions. `0` independent, `1`
  unstructured.

- NAA_region_corr_pars:

  Array `[pop, n_regions(n_regions-1)/2, sex]` of unconstrained
  parameters for the region correlation.

- NAA_re_pop, NAA_re_sex:

  Integer codes for the structure across populations and across sexes.
  `0` independent, `1` unstructured.

- NAA_pop_corr_pars, NAA_sex_corr_pars:

  Numeric vectors of unconstrained parameters for those correlations,
  one per pair. Both are global to the model rather than varying over
  the other margins, which is what keeps a two-level margin at exactly
  one parameter.

## Value

Scalar negative log likelihood.

## Details

The state is a level rather than a deviation, so the prediction is
subtracted here instead of multiplying a deviation onto a value the
dynamics still computes. That is what makes the random-effects Hessian
block-tridiagonal in year: every term is supported on a two-year block,
because \\\eta_y\\ depends on the states in years \\y\\ and \\y-1\\ and
on nothing earlier.

Only the independent form admits a standard deviation that varies cell
by cell. Every other structure is separable or Markov in a margin, and a
per-cell variance is neither, so the setup function holds the year and
age standard deviation blocks to one apiece whenever a correlated form
is chosen and this function reads one standard deviation per population,
region and sex.
