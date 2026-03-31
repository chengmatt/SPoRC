# Title Get Movement Process Error Likelihoods

Title Get Movement Process Error Likelihoods

## Usage

``` r
Get_move_PE_loglik(
  PE_model,
  PE_pars,
  move_devs,
  map_move_devs,
  do_recruits_move,
  adjacency_collapsed,
  move_type
)
```

## Arguments

- PE_model:

  Process error model values

- PE_pars:

  Process error parameters

- move_devs:

  Deviations

- map_move_devs:

  movement deviations to share

- do_recruits_move:

  Whether recruits move (0, don't move, 1 move)

- adjacency_collapsed:

  Adjacency matrix collapsed w/o retention

- move_type:

  Movement type (0 == unstructed; all regions connected, 1 == ctmc)

## Value

numeric value of log likelihood (in positive space)
