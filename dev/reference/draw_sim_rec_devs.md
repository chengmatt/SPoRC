# Draw the recruitment deviations for a simulation year

Writes each population and region's log recruitment deviation into
`sim_env$ln_RecDevs[, , y, sim]`. Deviation sharing follows
[`generate_initial_age_structure`](https://chengmatt.github.io/SPoRC/dev/reference/generate_initial_age_structure.md):
one draw per population when `n_pop > 1`, or one per region when
`n_pop = 1` under local density dependence. Populations with `R0 = 0`
get zero deviations, and `sigma_idx` picks the natal region's
`ln_sigmaR` for the bias correction. `RecDevs_model` sets what the draw
is centered on: zero for independent deviations, the previous year's for
a random walk, and `RecDevs_rho` times it for an AR1. Only the
independent draws are bias corrected, a walk's deviation not being mean
zero. A year covered by `Rec_input` draws nothing, and one covered by
`ln_RecDevs_input` reads its deviations from there.

## Usage

``` r
draw_sim_rec_devs(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place: `$ln_RecDevs`.

## Value

`invisible(NULL)`; everything is modified by reference within `sim_env`.

## Details

Under `rec_lag = 0`,
[`run_annual_cycle`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md)
draws year `y`'s deviations at the end of year `y - 1`, so catch advice
for year `y` can be converted to F knowing its recruits.
