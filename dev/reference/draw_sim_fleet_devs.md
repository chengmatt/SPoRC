# Draw fleet deviations for every replicate and rebuild the fleet arrays

F deviations under `Get_Fdev_PE_loglik`'s process, discard mortality
deviations iid, and selectivity deviations under `Get_PE_loglik`'s, each
kept at the fit over `n_cond_yrs` and drawn in the fitted years after
them. F and discard mortality are rebuilt in the cells whose deviation
was drawn, as `exp(ln_F_mean + ln_F_devs)` and
`plogis(logit_dmr_mean + logit_dmr_devs)`, and selectivity through
`Get_Selex_Array` over every fitted year. The draws are kept as
`ln_F_devs`, `logit_dmr_devs` and `ln_<type>sel_devs` with
`ln_<type>sel_bin_devs`, the replicate dim last.

## Usage

``` r
draw_sim_fleet_devs(sim_env)
```

## Arguments

- sim_env:

  Simulation environment holding what
  [`Setup_Sim_Fleet_Devs`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fleet_Devs.md)
  stored, `n_yrs`, `n_sims` and `n_cond_yrs` (`NULL` or `0` draws every
  fitted year).

## Value

`invisible(NULL)`; `sim_env` is modified in place.

## Details

Years past the fit, a closed loop's projection, draw selectivity
deviations on from the last fitted year under the same process
([`extend_devs_past_fit`](https://chengmatt.github.io/SPoRC/dev/reference/extend_devs_past_fit.md)),
and the projection years' selectivity is rebuilt over the longer span
and scaled so the fitted years keep the fit's standardization. Those
draws are kept as `ln_<type>sel_devs_proj`. A bicubic surface has no
year weights past the fit and keeps its last fitted year's selectivity.
