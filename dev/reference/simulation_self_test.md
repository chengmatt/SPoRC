# Run a simulation self-test of a fitted RTMB estimation model

Validates model performance by: (1) generating `n_sims` new datasets
from the fitted model parameters using
[`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
(2) re-fitting the estimation model to each simulated dataset, and (3)
storing user-specified report quantities for comparison against true
values. Supports sequential or parallel execution via
`future`/`future.apply`. Likelihood weights from the original fit are
propagated into the simulation (e.g., ISS scaled by `Wt_FishAgeComps`;
`ObsSrvIdx_SE` divided by `sqrt(Wt_SrvIdx)`), and the data weights are
reset to 1 when refitting. A penalty weight (`Wt_Rec`, `Wt_Init_Rec`,
`Wt_F`, `Wt_D`, `*_pe_wt`) is kept, the operating model drawing that
process at its sd over the root of the weight, so each refit is the
fit's own estimator. A weighted penalty whose sd is estimated comes back
at the sd over the root of the weight, a weighted penalty being a
density only up to its normalizing constant. Failed replicates are
silently stored as `NA`. The operating model takes the fit's
`bias_correct_pe`, `bias_correct_oe` and `sigmaR_switch`, so both sides
center process and observation error alike.

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
  obj = NULL,
  n_cores = NULL,
  output_path = NULL,
  what = c("SSB", "Rec"),
  what_par = NULL,
  perfect_data = FALSE,
  sim_type = c("conditional", "joint"),
  prior_means = c("assessment", "truth", "off")
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

  `sdreport` object from the fitted model, used to extract optimized
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

- obj:

  Fitted object the self test is run from. Needed under
  `sim_type = "joint"`, which draws at its parameter vector and reports
  through it. Default `NULL`.

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

- what_par:

  Character vector. Names of parameters (keys of `parameters`) to
  extract and store from each replicate, read off the refit's own
  parameter list so that mapped elements come back at the values the map
  gave them. An error is raised if any name is not found in
  `parameters`. Default `NULL`, which stores none.

- perfect_data:

  Logical. Whether to shrink the observation error before simulating,
  sds to 0.001 and sample sizes to 1e6, leaving process error alone. A
  correct model whose data identify the population then returns the
  operating model to several decimals. A state-space model need not:
  exact catch and survey indices at age still leave each age's scale to
  cohort continuity under process error, so recovery there is to a few
  percent. A process error sd conditioned on the fit comes back low,
  since the fitted value includes the posterior variance of its own
  deviations and data this clean remove it; on NEA cod the numbers at
  age sd returns 0.12 against 0.19. The truth for an estimated
  observation sd is the 0.001 the data were drawn at. Default `FALSE`.

- sim_type:

  Character. Which self test to run. `"conditional"` (default) runs
  every replicate at the fitted parameters and every process deviation
  at the fit's estimate, so one population covers every replicate and
  only the observations are new. It measures how well the data determine
  this history; a process sd comes back low, since the fit's deviations
  are shrunk estimates that vary less than the sd describes. `"joint"`
  draws each replicate's parameters from the fit's uncertainty
  (`sd_rep$jointPrecision`) and then every process fresh at them:
  recruitment and the initial ages from their penalty
  ([`rec_devs_past_fit`](https://chengmatt.github.io/SPoRC/dev/reference/rec_devs_past_fit.md),
  [`init_devs_past_fit`](https://chengmatt.github.io/SPoRC/dev/reference/init_devs_past_fit.md)),
  the numbers at age, selectivity, catchability, movement, growth and a
  linked dsem from their processes, and F and discard mortality
  deviations when the fit integrates them out (fixed effect F deviations
  come from the parameter draw). Each replicate is a new population with
  its own truth, so the self test measures bias in every estimate,
  process sds included, and with `do_sdrep = TRUE` whether the standard
  errors cover the truth. A fit with no random effects has no joint
  precision, and the fixed effect covariance is inverted in its place.

  Joint moves F, both selectivities, catchability, natural mortality,
  weight and size at age, movement, steepness, sex ratio, recruitment,
  the initial deviations, the numbers at age and a linked dsem.
  Observation error and the composition parameters stay at the fit,
  having no replicate dim.

- prior_means:

  Where each refit's priors are centered (see
  [`self_test_priors`](https://chengmatt.github.io/SPoRC/dev/reference/self_test_priors.md)).
  `"assessment"` (default) keeps them as the assessment gives them, so
  every replicate is pulled toward a prior the fit sits away from;
  `"truth"` centers each catchability, natural mortality, steepness, R0
  and selectivity prior on the value the replicate ran on, keeping its
  sd, and puts the mode of each Dirichlet (movement, recruitment
  apportionment) and mean and sd beta (stray rate, tag reporting) there,
  keeping its concentration; `"off"` turns every prior off.

## Value

Named list with one element per entry in `what` and then one per entry
in `what_par`, each an array with the last dimension indexing simulation
replicates (via `simplify2array`). If `do_sdrep = TRUE`, an additional
element `"sd_rep"` contains a list of `sdreport` objects (or `NA` for
failed replicates). A final element `"truth"` holds each replicate's
true values for the same names, the operating model's own wherever it
holds the quantity
([`om_truth`](https://chengmatt.github.io/SPoRC/dev/reference/om_truth.md)),
so replicates drawn under `sim_type = "joint"` are each compared with
what they ran on. Observation error, composition, at-age correlation and
tag loss parameters take the fit's values, which every replicate runs
at.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Fleet_Devs()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fleet_Devs.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
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

# parameters drawn from the fit's uncertainty, each replicate compared with its own truth
sd_rep <- RTMB::sdreport(fit, getJointPrecision = TRUE)
res <- simulation_self_test(
  data = fit$data, parameters = par, mapping = map, random = NULL,
  rep = fit$rep, sd_rep = sd_rep, obj = fit, n_sims = 100,
  what = "SSB", what_par = "ln_global_R0", sim_type = "joint"
)
rel_err <- (res$SSB - res$truth$SSB) / res$truth$SSB
} # }
```
