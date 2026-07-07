# Compute Selectivity Process Error Log-Likelihood (Positive Scale)

Calculates the positive log-likelihood contribution for selectivity
process error deviations under a variety of temporal/spatiotemporal
structures.

## Usage

``` r
Get_sel_PE_loglik(
  PE_model,
  PE_pars,
  ln_devs,
  map_sel_devs,
  min_sel_devs_shared_bins
)
```

## Arguments

- PE_model:

  Integer specifying the process error structure:

  - 1 = IID: deviations drawn independently as \\N(0, \sigma^2)\\.

  - 2 = Random walk: deviations follow a first-order random walk
    initialized with a diffuse prior (\\\sigma = 5\\) at `y = 1`.

  - 3 = 3D GMRF with marginal variance parameterization.

  - 4 = 3D GMRF with conditional variance parameterization.

  - 5 = Separable 2D AR(1) across bins and years.

- PE_pars:

  Array of process error parameters dimensioned
  `[1, par_index, sex, 1]`. The `par_index` slot meaning depends on
  `PE_model`:

  - Models 1–2: `[1,1,s,1]` = log standard deviation (\\\log \sigma\\)
    for sex `s`, indexed by bin/age.

  - Models 3–4: `[1,1,s,1]` = unconstrained partial correlation by
    age/bin; `[1,2,s,1]` = unconstrained partial correlation by year;
    `[1,3,s,1]` = unconstrained partial correlation by cohort;
    `[1,4,s,1]` = log variance.

  - Model 5: `[1,1,s,1]` = unconstrained bin correlation (transformed
    via \\2/(1+e^{-2x})-1\\); `[1,2,s,1]` = unconstrained year
    correlation; `[1,4,s,1]` = log standard deviation.

- ln_devs:

  Array of log-scale selectivity deviations dimensioned
  `[1, year, bin, sex, 1]`.

- map_sel_devs:

  Integer array dimensioned `[fleet, year, bin, sex]` mapping deviations
  to unique estimated parameters. Shared deviations carry the same
  integer value; `NA` entries are treated as fixed and excluded from
  likelihood evaluation.

- min_sel_devs_shared_bins:

  Integer vector. Indices of the reference (minimum) bin within each
  shared deviation group, used to subset the bin dimension when
  evaluating GMRF or 2D AR(1) likelihoods (PE models 3-5). When no bin
  sharing is specified, defaults to `1:n_bins` (i.e., all bins are
  included).

## Value

Numeric scalar: the positive log-likelihood contribution from
selectivity process error. Negated externally to form the negative
log-likelihood.

## Details

The function supports:

- IID process error

- Random walk process error

- 3D Gaussian Markov Random Field (GMRF) models (marginal or conditional
  variance)

- Separable 2D AR(1) models

**Note:** The returned value is on the *positive* log-likelihood scale.
It must be negated to obtain a negative log-likelihood contribution,
which is handled outside this function.
