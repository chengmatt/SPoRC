# Draw one replicate's catchability deviations and apply them

Called at the start of each replicate. Conditioning years keep the fit's
own deviations, a cell a dsem drew keeps that value, and every other
year is this replicate's own innovation, so a random walk or ar1
continues from the fit's last deviation rather than restarting at zero.

## Usage

``` r
draw_sim_q_devs(sim, sim_env)
```

## Arguments

- sim:

  Replicate index.

- sim_env:

  Simulation environment.

## Value

Nothing. `sim_env$fish_q` and `sim_env$srv_q` are scaled by the
deviations and the deviations themselves are stored.

## Details

A fleet a dsem wrote for is drawn whatever its own `<prefix>_q_model`
says, so the series reaches catchability rather than an array nothing
reads.
