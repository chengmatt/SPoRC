# Draw movement deviations for every replicate and rebuild movement

The first `n_cond_yrs` years of every replicate take the fit's own
deviations, or the replicate's own draw of them under a joint self test,
and the years after them are drawn from the process the estimation model
penalizes, at that replicate's parameters: a stationary AR1 over ages
where that dim is `"ar1"`, an unstructured correlation across
populations, seasons or sexes where a dim is `"us"`, independent
otherwise, with the conditional sd `exp(move_pe_pars[..., 1])` raised to
the marginal by `1 / sqrt(1 - rho^2)` for each `"ar1"` dim, as the
estimation model reads it. An AR1 over years continues from the last
conditioned year rather than restarting, a shared year deviation
(`"none"`) continues unchanged, and with no conditioned years the whole
series is drawn from its stationary distribution. One value is drawn per
block of each dim and written into every level of that block, as the
estimation model's map ties them; inactive levels stay at zero,
projection years are always active, and a cell the map holds as `NA`
stays at zero. Movement is then rebuilt through `derive_sim_movement`.

## Usage

``` r
draw_sim_move_devs(sim_env)
```

## Arguments

- sim_env:

  Simulation environment holding what `Setup_Sim_Movement` stored,
  `n_yrs`, `n_sims` and `n_cond_yrs` (`NULL` or `0` draws every year).

## Value

`invisible(NULL)`; `sim_env` is modified in place.
