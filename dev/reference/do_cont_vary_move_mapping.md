# Map continuous movement deviation and process-error parameters

Internal helper called by
[`Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md)
to construct the TMB/RTMB factor maps for `move_devs` (iid deviations on
the movement logit or log-rate surface) and `move_pe_pars`
(process-error variance parameters). Deviations are only activated when
the model is spatial (`n_regions > 1`), continuous variation is
requested (`cont_vary_movement > 0`), and movement is estimated
(`use_fixed_movement == 0`). For CTMC movement, deviations are only
assigned to region pairs that are connected in the adjacency matrix. The
resulting integer map is also stored as `$data$map_move_devs` for use in
the C++ template.

## Usage

``` r
do_cont_vary_move_mapping(
  input_list,
  cont_vary_movement,
  Movement_cont_pe_pars_spec
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- cont_vary_movement:

  Character string specifying the deviation structure. One of `"none"`,
  `"iid_y"`, `"iid_a"`, `"iid_y_a"`, `"iid_y_a_s"`, `"iid_y_seas_a_s"`,
  or the population-specific analogues `"iid_p_y"`, `"iid_p_a"`,
  `"iid_p_y_a"`, `"iid_p_y_a_s"`, `"iid_p_y_seas_a_s"`. Dimensions
  present in the string receive unique estimation indices; absent
  dimensions share a single index. `"none"` maps all deviations to `NA`.

- Movement_cont_pe_pars_spec:

  Character string controlling estimation of the process-error variance
  for movement deviations. One of:

  `"none"` or `"fix"`

  :   All `move_pe_pars` held fixed (mapped to `NA` or at starting
      values).

  `"est_shared"`

  :   Single variance parameter shared across all dimensions (all
      elements mapped to index 1).

  `"est_all"`

  :   All `move_pe_pars` estimated independently with dimensions
      `[n_pop × n_regions × n_seas × n_ages × n_sexes]`.

## Value

The input `input_list` with three entries updated:

- `$map$move_devs`:

  Factor vector for movement deviations. Active cells receive sequential
  integer indices; non-adjacent CTMC pairs and inactive configurations
  are `NA`.

- `$data$map_move_devs`:

  Integer array (same dimensions as `$par$move_devs`) storing the
  numeric version of the factor map for use in the C++ objective
  function.

- `$map$move_pe_pars`:

  Factor vector for process-error variance parameters, following
  `Movement_cont_pe_pars_spec`.
