# Run the annual cycle for a single simulation year

Orchestrates the complete annual sequence of operating model processes
for year `y` and simulation replicate `sim`: initialises age structure
and generates first-year recruitment at `y = 1`; applies population
dynamics (movement, mortality, biomass); generates fishery catches,
indices, and compositions; generates survey indices and compositions;
releases conventional tags; generates fishery tag recaptures (when any
`use_conv_fish_tagging = 1`); and generates recruitment for the
following year (`y + 1`) when `y < n_yrs`.

## Usage

``` r
run_annual_cycle(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment created by
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
  and passed by reference. All annual-cycle helper functions modify this
  environment in place.

## Value

`invisible(NULL)`.

## Details

The two standalone
[`generate_recruitment()`](https://chengmatt.github.io/SPoRC/dev/reference/generate_recruitment.md)
calls described above (at `y = 1` and for `y + 1`) only run when
`rec_lag != 0`. For `rec_lag = 0` (age-0 recruitment), recruitment for
year `y` depends on year `y`'s own SSB, which isn't known until
[`apply_pop_dy`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md)
reaches `spawn_seas` within that year -
[`generate_recruitment()`](https://chengmatt.github.io/SPoRC/dev/reference/generate_recruitment.md)
is called from inside
[`apply_pop_dy()`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md)
instead, once that SSB is available.

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
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
