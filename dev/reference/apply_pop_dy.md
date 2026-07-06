# Apply population dynamics within a simulation year

Executes the full within-year population dynamics loop for year `y`:
seasonal recruitment apportionment (seasons 2+), movement, Baranov
catch-equation mortality, age advancement into the following year, and
spawning-season biomass calculations (total biomass, SSB, dynamic
\\B_0\\, and effective SSB for multi-population natal homing). Both
fished (`NAA`) and unfished (`NAA0`) trajectories are tracked in
parallel. For single-season multi-population models,
`sgl_seas_spawning_movement` is applied to `NAA` and `NAA0` prior to
computing spawning biomass quantities. Single-sex models have SSB and
\\B_0\\ multiplied by 0.5 to obtain female-only spawning biomass.

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

  Simulation environment created by
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).
  Modified in place: `$ZAA`, `$NAA`, `$NAA0`, `$NAA_bef`, `$NAA_aft`,
  `$Total_Biom`, `$SSB`, `$Dynamic_SSB0`, and `$eff_SSB` are updated.

## Value

`invisible(NULL)`. All modifications are made by reference within
`sim_env`.

## Details

Pre- and post-movement snapshots are stored in `NAA_bef` and `NAA_aft`
respectively. Movement is only applied when `n_regions > 1`; recruits
(`a = 1`) are excluded from movement when `do_recruits_move = 0`.

When `rec_lag == 0` (age-0 recruitment), this year's recruitment can't
be known until `spawn_seas` is reached (it depends on this year's own
SSB), so
[`generate_recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/generate_recruitment.md)
is called from inside this function at `seas == spawn_seas` instead of
beforehand - see the "rec_lag == 0" block below, which mirrors the
equivalent restructuring in the estimation model (`SPoRC_rtmb.R`).
