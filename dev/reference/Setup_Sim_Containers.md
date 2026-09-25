# Initialize output containers for the operating model simulation

Allocates the zero-filled arrays the simulation writes into, sized off
the dimensions in `sim_list`. Call after
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md)
and before any operating model dynamics.

## Usage

``` r
Setup_Sim_Containers(sim_list)
```

## Arguments

- sim_list:

  A simulation list returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
  whose `n_pop`, `n_regions`, `n_yrs`, `n_seas`, `n_ages`, `n_sexes`,
  `n_sims`, `n_fish_fleets`, `n_srv_fleets`, `n_obs_ages` and `n_lens`
  size every container.

## Value

`sim_list` with the containers added.

Biological: `$NAA`
`[n_pop × n_regions × (n_yrs+1) × n_seas × n_ages × n_sexes × n_sims]`,
whose extra year holds the initial conditions and advances the
population through the final year; `$NAA_bef` and `$NAA_aft`, the
numbers before and after fishing mortality; `$NAA0`, the unfished
numbers, for dynamic \\B_0\\; `$ZAA`, total mortality at age over
`n_yrs`; and `$Rec`, `$SSB`, `$Dynamic_SSB0` and `$Total_Biom`
`[n_pop × n_regions × n_yrs × n_sims]`, with `$eff_SSB`
`[n_pop × n_yrs × n_sims]`, `$ln_RecDevs` on the SSB dims and
`$ln_InitDevs` `[n_pop × n_regions × (n_ages - 1) × n_sexes × n_sims]`.

Fishery: `$ObsCatch` and `$TrueCatch`
`[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]`, with
`$ObsFishIdx`, `$TrueFishIdx`, `$ObsDiscard` and `$TrueDiscard` on the
same dims; `$ObsFishAgeComps` and `$ObsFishAgeComps_discard` with
`n_obs_ages × n_sexes` before the fleet dim, and `$ObsFishLenComps` and
`$ObsFishLenComps_discard` with `n_lens × n_sexes`. Each has a `_pop`
counterpart with a leading `n_pop`. The true catch and discards at age
and length are `$CAA`, `$DAA`
`[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims]`
and `$CAL`, `$DAL` with `n_lens` in place of `n_ages`.

Survey: `$ObsSrvIdx` and `$TrueSrvIdx`
`[n_regions × n_yrs × n_seas × n_srv_fleets × n_sims]`,
`$ObsSrvAgeComps` and `$ObsSrvLenComps` with the bin and sex dims before
the fleet dim, each with a `_pop` counterpart, and the true `$SrvIAA`
`[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]`
and `$SrvIAL` with `n_lens` in place of `n_ages`.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
