# Apply population dynamics within a simulation year

Runs the within-year loop for year `y`: seasonal recruitment
apportionment from season two on, movement, Baranov mortality, age
advancement into the following year, and the spawning-season biomass
quantities (total biomass, SSB, dynamic \\B_0\\ and effective SSB under
natal homing). The fished and unfished trajectories are tracked
together, with the snapshots before and after movement stored in
`NAA_bef` and `NAA_aft`. Movement runs only when `n_regions > 1`, and
recruits are left out of it when `do_recruits_move = 0`. With one season
and several populations, `sgl_seas_spawning_movement` is applied to both
trajectories before the biomass quantities; with one sex, SSB and
\\B_0\\ are halved to give female-only spawning biomass.

## Usage

``` r
apply_pop_dy(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place: `$ZAA`, `$NAA`, `$NAA0`, `$NAA_bef`, `$NAA_aft`,
  `$Total_Biom`, `$SSB`, `$Dynamic_SSB0` and `$eff_SSB`.

## Value

`invisible(NULL)`; everything is modified by reference within `sim_env`.

## Details

Under `rec_lag == 0` this year's recruitment is not knowable until
`spawn_seas`, since it depends on this year's own SSB, so
[`generate_recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/generate_recruitment.md)
is called from inside this function at `seas == spawn_seas`, mirroring
the estimation model.
