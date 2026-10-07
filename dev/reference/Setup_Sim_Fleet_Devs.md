# Fishing mortality, discard mortality and selectivity deviations in the operating model

Stores what each replicate needs to draw these deviations and rebuild
the fleet arrays from them, and
[`draw_sim_fleet_devs`](https://chengmatt.github.io/SPoRC/dev/reference/draw_sim_fleet_devs.md)
does so in `Setup_sim_env`, keeping the fit's deviations over the
simulation list's `n_cond_yrs` and drawing the fitted years after them.

## Usage

``` r
Setup_Sim_Fleet_Devs(
  sim_list,
  data,
  pars,
  random = NULL,
  pars_by_sim = NULL,
  F_devs_draw = c("random", "all")
)
```

## Arguments

- sim_list:

  Simulation list with `n_yrs`, `n_sims`, `Fmort`, `dmr` and the
  selectivity arrays already set.

- data:

  Data list of the fit.

- pars:

  Parameter list at the fitted values.

- random:

  Character vector of the fit's random effects. Default `NULL`.

- pars_by_sim:

  Optional list of one parameter list per replicate, for replicates that
  each run on their own parameter draw (`sim_type = "joint"` in
  [`simulation_self_test`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)).
  `NULL` (default) gives every replicate `pars`.

- F_devs_draw:

  `"random"` (default) draws F and discard mortality deviations only
  when `random` integrates them out; `"all"` also draws them when the
  fit holds them as fixed effects, from their penalty.

## Value

`sim_list` with `fleet_devs_on` (which of `F`, `dmr`, `fish`, `ret` and
`srv` are drawn), `fleet_dev_data` and `fleet_dev_pars` added, or
unchanged when nothing is drawn.

## Details

Selectivity deviations are drawn wherever the fit varies selectivity,
under the process
[`Get_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_PE_loglik.md)
penalizes, as growth and movement deviations are. F and discard
mortality deviations are drawn by default only when the fit integrates
them out (`random`): as fixed effects their penalty is usually a loose
constraint, `sigmaF = 1` by default, and drawing from it would replace
the fitted F with noise. `F_devs_draw = "all"` draws them as fixed
effects too, from the same penalty. Either way they are kept at the fit
when the penalty is off or centered on the deviations' own mean, neither
of which is a density to draw from.

A penalty weighted by `Wt_F`, `Wt_D` or `*_pe_wt` is the density of the
same process with its variance divided by the weight, and that is what
is drawn; a self test keeps the weight in its refits.

Selectivity at length (`*_selex_type = 1`) is drawn and rebuilt at
length, and its selectivity at age is the curve read through each
replicate's own size-age transition matrix, by the growth rebuild when
the operating model rebuilds growth and directly otherwise. A penalty
weight `w` (`*_pe_wt`) scales the penalty to the density of a process
whose variance is divided by `w`, which is what is drawn. When the
operating model runs past the fit, as a closed loop does, selectivity
deviations are drawn in those years too, on from the fitted ones (see
[`draw_sim_fleet_devs`](https://chengmatt.github.io/SPoRC/dev/reference/draw_sim_fleet_devs.md)).
F and discard mortality never are, since a closed loop sets F in its
projection years by its control rule.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
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
