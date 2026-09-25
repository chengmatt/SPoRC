# Initialize age structure for a simulation replicate

Draws or reads the initial age deviations and calls
[`Get_Init_NAA`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Init_NAA.md)
for the fished and unfished equilibrium numbers at age in year 1, season
1, writing them into `NAA` and `NAA0`. Called once per replicate at
`y = 1` by
[`run_annual_cycle`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md).

## Usage

``` r
generate_initial_age_structure(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index, which must be `1`.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place: `$ln_InitDevs`, `$NAA[,,1,1,,,sim]` and
  `$NAA0[,,1,1,,,sim]`.

## Value

`invisible(NULL)`; everything is modified by reference within `sim_env`.

## Details

Sharing follows the estimation model: one draw per population when
`n_pop > 1`, or one per region when `n_pop = 1` and `init_dd = 0`.
Across sexes it follows `InitDevs_sex_spec`, with `"est_shared_s"`
(default) drawing one curve for every sex and `"est_all"` drawing each
its own. An `ln_InitDevs_input` in the environment is used directly
rather than drawn. Populations with `R0 = 0` get zero deviations, and
the equilibrium solver runs `n_ages × 5` iterations.
