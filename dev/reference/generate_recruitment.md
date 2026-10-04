# Generate recruitment for a simulation year

Takes deterministic recruitment from `sim_det_recruitment`, multiplies
it by the year's lognormal deviations from `draw_sim_rec_devs`,
apportions it across sexes and seasons, and writes it into the age-one
slot of `sim_env$NAA`, with `NAA0` synchronized to match. A `Rec_input`
covering year `y` overrides the draw entirely.

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
