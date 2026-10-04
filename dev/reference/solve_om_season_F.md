# Solve One Season's F For A Block Of Regions

Tapes log catch in the cells being solved as a function of their log F
with [`RTMB::MakeTape`](https://rdrr.io/pkg/RTMB/man/Tape.html) and
solves them jointly with `solve_log_catch_newton`, since catch in a cell
falls as other fleets in the same region fish harder.

## Usage

``` r
solve_om_season_F(
  F_seas,
  free,
  target_seas,
  N,
  y,
  seas,
  sim,
  sim_env,
  target_units,
  spawn_recruits,
  catch_f_max,
  catch_tol,
  catch_max_iter
)
```

## Arguments

- F_seas:

  Numeric matrix `[n_regions, n_fish_fleets]`. This season's F so far;
  cells outside `free` are kept fixed.

- free:

  Integer vector. Cells to solve, counted down regions first.

- target_seas:

  Numeric matrix `[n_regions, n_fish_fleets]`.

- N, y, seas, sim, sim_env, target_units, spawn_recruits:

  Passed to `om_season_catch`.

- catch_f_max, catch_tol, catch_max_iter:

  Solver settings.

## Value

Named list with `F_seas`, solved in the `free` cells, and `capped`, the
cells left at `catch_f_max` because the target is not reachable there.
