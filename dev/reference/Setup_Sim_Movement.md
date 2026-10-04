# Movement deviations in the operating model

Stores what each replicate needs to rebuild movement from its own
deviations: the fit's switches, deviation map, surfaces, active years,
seasons and ages, process error and correlation parameters, and the
`Get_Movement` arguments, and the fit's own deviations. With the fit's
own penalty on the deviations,
[`draw_sim_move_devs`](https://chengmatt.github.io/SPoRC/dev/reference/draw_sim_move_devs.md)
keeps the fit's values over the simulation list's `n_cond_yrs` and draws
the years after them in `Setup_sim_env`; with a dsem holding or linking
them,
[`Setup_Sim_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_DSEM.md)
draws them and the arguments stored here rebuild movement. Nothing is
stored when the fit has no movement deviations or movement is fixed.

## Usage

``` r
Setup_Sim_Movement(sim_list, data, pars, pars_by_sim = NULL)
```

## Arguments

- sim_list:

  Simulation list with `n_yrs`, `n_sims` and `expm_nsub`.

- data:

  Data list of the fit.

- pars:

  Parameter list at the fitted values.

- pars_by_sim:

  Optional list of one parameter list per replicate, for replicates that
  each run on their own parameter draw (`sim_type = "joint"` in
  [`simulation_self_test`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)).
  Each replicate's movement is then rebuilt from its own parameters,
  deviations and process error. `NULL` (default) gives every replicate
  `pars`.

## Value

`sim_list` with `move_args` and the fields above added.
