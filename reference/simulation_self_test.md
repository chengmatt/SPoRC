# Conduct a Simulation Self Test

This function runs a self test of the fitted RTMB model by simulating
new datasets under the fitted parameters, refitting the model, and
comparing estimated outputs to the true values used for simulation. It
can be run sequentially or in parallel.

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

  A list containing model data from an RTMB object.

- parameters:

  A list of fitted parameter values from an RTMB object.

- mapping:

  A list specifying parameter mappings from an RTMB object.

- random:

  Character vector specifying random effects.

- rep:

  A list of report values from an RTMB object (\`\$rep\`).

- sd_rep:

  An \`sdreport\` object from RTMB summarizing parameter uncertainty.

- n_sims:

  Integer. Number of simulation replicates to run.

- newton_loops:

  Integer. Number of Newton loops used in model fitting (default:
  \`3\`).

- do_sdrep:

  Logical. If \`TRUE\`, compute \`sdreport\` for each fitted replicate
  (default: \`FALSE\`).

- do_par:

  Logical. If \`TRUE\`, run simulations in parallel (default:
  \`FALSE\`).

- n_cores:

  Integer. Number of cores to use for parallelization (default: \`NULL\`
  = detect automatically).

- output_path:

  Optional file path. If provided, the simulated datasets are written to
  this location.

- what:

  Character vector. Names of report elements in \`rep\` to extract and
  store for each replicate.

## Value

A list with elements corresponding to the requested \`what\` values,
each containing an array of simulation results across replicates. If
\`do_sdrep = TRUE\`, an additional element \`"sd_rep"\` is included with
the list of \`sdreport\` objects (or \`NA\` if a replicate fails).

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Run a simple self test with 10 simulations, extracting SSB
res <- simulation_self_test(
  data = model$data,
  parameters = model$parameters,
  mapping = model$mapping,
  random = model$random,
  rep = model$rep,
  sd_rep = model$sd_rep,
  n_sims = 10,
  what = "SSB"
)

str(res$SSB) # look at simulated SSB arrays
} # }
```
