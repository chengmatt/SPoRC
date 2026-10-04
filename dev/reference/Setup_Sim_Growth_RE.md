# Growth deviations in the operating model

Stores what each replicate needs to rebuild growth from its own
deviations: the fit's time variation codes, the semi-parametric form and
its ages, the deviation maps, the process error parameters, the
`Get_Growth` arguments and, under selectivity at length, the report's
selectivity at length the rebuilt keys are read through. With the fit's
own penalty on the deviations,
[`draw_sim_growth_devs`](https://chengmatt.github.io/SPoRC/dev/reference/draw_sim_growth_devs.md)
keeps the fit's values over the simulation list's `n_cond_yrs` and draws
the years after them in `Setup_sim_env`. With a dsem holding or linking
any growth deviation,
[`Setup_Sim_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_DSEM.md)
builds both arrays itself and nothing is stored here, so the deviations
of a parameter outside the arrows keep the fit's values. Nothing is
stored either when growth is data or has no deviations.

## Usage

``` r
Setup_Sim_Growth_RE(
  sim_list,
  data,
  pars,
  rep = NULL,
  pars_by_sim = NULL,
  rep_by_sim = NULL
)
```

## Arguments

- sim_list:

  Simulation list with `n_yrs` and `n_sims`.

- data:

  Data list of the fit.

- pars:

  Parameter list at the fitted values.

- rep:

  Report of the fit, needed under selectivity at length.

- pars_by_sim, rep_by_sim:

  Optional lists of one parameter list and one report per replicate, for
  replicates that each run on their own parameter draw
  (`sim_type = "joint"` in
  [`simulation_self_test`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)).
  Each replicate's growth is then rebuilt from its own parameters and
  deviations. `NULL` (default) gives every replicate `pars` and `rep`.

## Value

`sim_list` with `growth_args` and the fields above added.
