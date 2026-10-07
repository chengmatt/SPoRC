# The priors a self test's refits read

A prior's mean is data the refit reads, but a self test does not redraw
it with the observations, so where the data pulled the fit away from a
prior every replicate is pulled back toward it, a bias of the design
rather than of the estimator. `"assessment"` keeps the priors as the
assessment gives them; `"truth"` moves the mean of every catchability,
natural mortality, steepness, R0 and selectivity prior to the value the
replicate ran on, keeping its sd, so the prior keeps its information and
loses its pull; `"off"` turns every prior off, which leaves a parameter
only a prior identifies unidentified. A Dirichlet (on movement and on
the recruitment apportionment) or a mean and sd beta (on the stray rate
and the tag reporting rate) is moved the same way, its mode put on the
value the replicate ran on and its concentration kept, so it has no
gradient there
([`dirichlet_mode_at`](https://chengmatt.github.io/SPoRC/dev/reference/dirichlet_mode_at.md),
[`beta_mode_at`](https://chengmatt.github.io/SPoRC/dev/reference/beta_mode_at.md));
one with no interior mode is left as given, as is the symmetric beta on
the tag reporting rate, which sits at one half by construction.

## Usage

``` r
self_test_priors(data, prior_means, pars, rep)
```

## Arguments

- data:

  Data list of the fit.

- prior_means:

  `"assessment"`, `"truth"` or `"off"`.

- pars, rep:

  The parameter list and report the replicate ran on.

## Value

`data` with its priors set.
