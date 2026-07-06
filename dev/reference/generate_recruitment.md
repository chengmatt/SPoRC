# Generate recruitment for a simulation year

Computes total recruitment for year `y` by first obtaining deterministic
expected recruitment from
[`Get_Det_Recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Det_Recruitment.md)
and then multiplying by lognormal deviations (bias-corrected).
Recruitment is apportioned across sexes and seasons and written into
`sim_env$NAA[p, r, y, 1, 1, s, sim]` (first season, age-1 slot).
Unfished NAA (`NAA0`) is synchronised to match fished NAA at
recruitment. If `Rec_input` exists in the environment and covers year
`y`, those values override the stochastic draw entirely.

## Usage

``` r
generate_recruitment(y, sim, sim_env, seas = 1)
```

## Arguments

- y:

  Integer. Year index for which recruitment is generated.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment created by
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).
  Modified in place: `$ln_RecDevs[p, r, y, sim]`, `$Rec[p, r, y, sim]`,
  `$NAA[p, r, y, seas, 1, s, sim]`, and
  `$NAA0[p, r, y, seas, 1, s, sim]` are updated.

- seas:

  Integer. Season this recruitment first enters the population in, using
  `rec_seas_prop[p, seas, sim]`. Default `1`, matching the classic
  `rec_lag >= 1` case where recruitment for the whole year is already
  known before season 1 starts. `rec_lag = 0` (age-0 recruitment)
  instead calls this with `seas = spawn_seas`, since that is the
  earliest season this year's own SSB (and hence recruitment) is
  knowable – see
  [`apply_pop_dy`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md).

## Value

`invisible(NULL)`. All modifications are made by reference within
`sim_env`.

## Details

Recruitment deviation sharing follows the same population/region logic
as
[`generate_initial_age_structure`](https://chengmatt.github.io/SPoRC/dev/reference/generate_initial_age_structure.md):
deviations are drawn once per population (`n_pop > 1`) or once per
region (`n_pop = 1`, local density-dependence). Populations with
`R0 = 0` receive zero deviations. `sigma_idx` selects the natal region's
`ln_sigmaR` for the bias-correction term.
