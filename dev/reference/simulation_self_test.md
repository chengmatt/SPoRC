# Run a simulation self-test of a fitted RTMB estimation model

Validates model performance by: (1) generating `n_sims` new datasets
from the fitted model parameters using
[`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
(2) re-fitting the estimation model to each simulated dataset, and (3)
storing user-specified report quantities for comparison against true
values. Supports sequential or parallel execution via
`future`/`future.apply`. Likelihood weights from the original fit are
propagated into the simulation (e.g., ISS scaled by `Wt_FishAgeComps`;
`ObsSrvIdx_SE` divided by `sqrt(Wt_SrvIdx)`); all weights are reset to 1
when re-fitting. Failed replicates are silently stored as `NA`.

## Usage

``` r
simulation_self_test(
  data,
  parameters,
  mapping,
  random,
  rep,
  sd_rep,
  n_sims,
  newton_loops = 3,
  do_sdrep = FALSE,
  do_par = FALSE,
  n_cores = NULL,
  output_path = NULL,
  what = c("SSB", "Rec")
)
```

## Arguments

- data:

  Named list of model data from a fitted RTMB object (`$data`).

- parameters:

  Named list of fitted parameter values (`$par` or equivalent).

- mapping:

  Named list of parameter factor maps (`$map`).

- random:

  Character vector of random effect names passed to RTMB.

- rep:

  Named list of report values from the fitted model (`obj$rep`).

- sd_rep:

  `sdreport` object from the fitted model, used to extract optimised
  parameter values in list format via `get_optim_param_list`.

- n_sims:

  Integer. Number of simulation replicates.

- newton_loops:

  Integer. Number of Newton refinement steps applied during re-fitting.
  Default `3`.

- do_sdrep:

  Logical. Whether to compute `sdreport` for each fitted replicate.
  Results stored as `$sd_rep` in the output list; failed `sdreport`
  calls stored as `NA`. Default `FALSE`.

- do_par:

  Logical. Whether to run replicates in parallel via
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html).
  Default `FALSE`.

- n_cores:

  Integer. Number of parallel workers. If `NULL` (default),
  `parallel::detectCores() - 1` is used.

- output_path:

  Character string. Path to save the simulated dataset RDS file. Passed
  to
  [`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md).
  Default `NULL`.

- what:

  Character vector. Names of report elements (keys of `rep`) to extract
  and store from each replicate. An error is raised if any name is not
  found in `rep`. Default `c("SSB", "Rec")`.

## Value

Named list with one element per entry in `what`, each an array with the
last dimension indexing simulation replicates (via `simplify2array`). If
`do_sdrep = TRUE`, an additional element `"sd_rep"` contains a list of
`sdreport` objects (or `NA` for failed replicates).

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- simulation_self_test(
  data = obj$data, parameters = obj$par, mapping = obj$map,
  random = obj$random, rep = obj$rep, sd_rep = obj$sd_rep,
  n_sims = 100, what = c("SSB", "Rec", "Fmort")
)
str(res$SSB)
} # }
```
