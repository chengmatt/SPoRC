# Refuse a catchability series the operating model cannot draw

A dsem series given an sd of zero is derived: it has no innovation of
its own, the objective works it out from its covariates, and the fitted
deviation parameter stays at zero.
[`draw_dsem_sim`](https://chengmatt.github.io/SPoRC/dev/reference/draw_dsem_sim.md)
draws the unknown cells from the field's precision, which a series with
no variance has no finite entry in, so every cell it has to draw comes
back `NaN`.

## Usage

``` r
check_q_dsem_drawable(sim_list)
```

## Arguments

- sim_list:

  Simulation list, or the environment
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
  built from one.

## Value

Nothing. Called for its refusal.

## Details

Catchability is the only process a derived series is allowed on, since
[`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md)
refuses one everywhere else, and the conditioning years are read from
the fit rather than drawn. So the one case with cells left to draw is a
derived catchability series on an operating model that runs past the
conditioning period, and that is what this refuses.
