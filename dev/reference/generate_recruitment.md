# Generate recruitment for a simulation year

Takes deterministic recruitment from
[`Get_Det_Recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Det_Recruitment.md),
multiplies it by lognormal deviations, apportions it across sexes and
seasons, and writes it into the age-one slot of `sim_env$NAA`, with
`NAA0` synchronized to match. A `Rec_input` covering year `y` overrides
the draw entirely.

## Usage

``` r
generate_recruitment(y, sim, sim_env, seas = 1)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place: `$ln_RecDevs`, `$Rec`, `$NAA` and `$NAA0`.

- seas:

  Integer. Season this recruitment first enters in, through
  `rec_seas_prop[p, seas, sim]`. Default `1`, the classic `rec_lag >= 1`
  case where the whole year's recruitment is known before season one.
  `rec_lag = 0` instead calls this with `seas = spawn_seas`, the
  earliest season this year's own SSB is knowable. See
  [`apply_pop_dy`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md).

## Value

`invisible(NULL)`; everything is modified by reference within `sim_env`.

## Details

Deviation sharing follows
[`generate_initial_age_structure`](https://chengmatt.github.io/SPoRC/dev/reference/generate_initial_age_structure.md):
one draw per population when `n_pop > 1`, or one per region when
`n_pop = 1` under local density dependence. Populations with `R0 = 0`
get zero deviations, and `sigma_idx` picks the natal region's
`ln_sigmaR` for the bias correction. `RecDevs_model` sets what the draw
is centered on: zero for independent deviations, the previous year's for
a random walk, and `RecDevs_rho` times it for an AR1. Only the
independent draws are bias corrected, a walk's deviation not being mean
zero.
