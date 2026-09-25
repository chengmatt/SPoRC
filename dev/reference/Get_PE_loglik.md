# Compute Process Error Log-Likelihood for a Deviation Surface (Positive Scale)

The positive log-likelihood of a surface of deviations indexed by year
and a second dim, under iid, random walk, 3D GMRF (marginal or
conditional variance) or separable 2D AR(1) process error. Selectivity
deviations use it over years and bins, growth's semi-parametric
deviations over years and ages, and a time-varying growth parameter over
years alone, a surface one column wide; the argument names read `bin`
for that second dim throughout. The caller negates the result.

## Usage

``` r
Get_PE_loglik(
  PE_model,
  PE_pars,
  ln_devs,
  map_sel_devs,
  map_sel_devs_full,
  min_sel_devs_shared_bins,
  rw_init_sigma = 5
)
```

## Arguments

- PE_model:

  Integer process error structure: `1` iid, `2` random walk with a
  diffuse prior at `y = 1`, `3` and `4` the 3D GMRF on the marginal or
  conditional variance, `5` the separable 2D AR(1) over bins and years.

- PE_pars:

  Array of process error parameters `[1, par_index, sex, 1]`, whose
  `par_index` slots depend on `PE_model`. Models 1 and 2 hold a log
  standard deviation in slot 1, indexed by bin. Models 3 and 4 hold the
  unconstrained partial correlations by bin, year and cohort in slots 1
  to 3 and a log variance in slot 4. Model 5 holds the unconstrained bin
  and year correlations in slots 1 and 2 and a log standard deviation in
  slot 4.

- ln_devs:

  Array of log-scale deviations `[1, year, bin, sex, 1]`.

- map_sel_devs:

  Integer array `[fleet, year, bin, sex]` mapping the deviations to
  estimated parameters. Shared deviations hold the same integer, and
  `NA` entries are fixed and left out of the likelihood.

- map_sel_devs_full:

  The same map across every unit this penalty is evaluated over, that
  unit being the first dim: regions for selectivity, populations by
  region for growth. A deviation shared over those units is one
  parameter appearing in each of their slices, and this function runs
  one unit at a time, so its contribution is divided by the number
  holding it. Without the split, a series shared over `n` units is
  penalized `n` times, an implicit \\\sigma / \sqrt{n}\\.

- min_sel_devs_shared_bins:

  Integer vector of the reference bin within each shared deviation
  group, used to subset the bin dim under process error models 3 to 5.
  Defaults to `1:n_bins` when no bin sharing is set.

- rw_init_sigma:

  Standard deviation given to the first year of a random walk. A number
  (5 by default) leaves that year effectively unconstrained; `NA` starts
  the walk at zero under its own sigma.

## Value

Numeric scalar, the positive log-likelihood, negated by the caller.
