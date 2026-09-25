# Draw dsem series for every replicate

Draws the year by series grid for all replicates up front, which is the
same in distribution as drawing year by year because the series do not
depend on the harvest, and keeps the same draws across management
procedures. A linked recruitment cell is drawn about minus half its
variance given the known cells
([`get_dsem_margvar`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_margvar.md),
with a moderating series at its mean) when `rec_bias_correct` is on.

## Usage

``` r
draw_dsem_sim(sim_env)
```

## Arguments

- sim_env:

  Simulation environment from `Setup_sim_env`, built from a list set up
  by `Setup_Sim_DSEM`.

## Value

`invisible(NULL)`; `sim_env` is modified in place.

## What it fills

`dsem_x_sim` `[year, series, sim]`, every linked cell of every linked
array from it with `dsem_drawn` marking those cells, and
`dsem_cov_obs_sim` `[year, covariate, sim]` (`NA` where not observed).
Growth and movement are then rebuilt from the arrays a drawn series went
into.
