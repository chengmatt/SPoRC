# Run the annual cycle for a single simulation year

Runs the operating model's processes for year `y` and replicate `sim` in
order: initial age structure and first-year recruitment at `y = 1`,
population dynamics, fishery catches, indices and compositions, survey
indices and compositions, tag releases, fishery tag recaptures when any
`use_conv_fish_tagging = 1`, and next year's recruitment when
`y < n_yrs`.

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

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  passed by reference and modified in place by every annual-cycle
  helper.

## Value

`invisible(NULL)`.

## Details

Those two standalone
[`generate_recruitment()`](https://chengmatt.github.io/SPoRC/dev/reference/generate_recruitment.md)
calls only run when `rec_lag != 0`. Under `rec_lag = 0` recruitment
depends on year `y`'s own SSB, which is not known until
[`apply_pop_dy`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md)
reaches `spawn_seas`, so it is called from inside
[`apply_pop_dy()`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md)
instead.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
