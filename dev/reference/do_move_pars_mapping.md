# Map unstructured Markov movement parameters

Internal helper called by
[`Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md)
to construct the TMB/RTMB factor maps for `move_pars` (unstructured
Markov transitions), `log_move_diffusion_pars` (CTMC diffusion), and
`move_preference_pars` (CTMC taxis). Under the unstructured Markov
formulation (`move_type = 0`), parameters within each block combination
share a common estimation index; for CTMC (`move_type = 1`), `move_pars`
is mapped entirely to `NA` and only the CTMC-specific arrays are
activated.

## Usage

``` r
do_move_pars_mapping(
  input_list,
  Movement_popblk_spec,
  Movement_ageblk_spec,
  Movement_yearblk_spec,
  Movement_sexblk_spec,
  Movement_seasblk_spec,
  use_fixed_movement
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions.

- Movement_popblk_spec:

  `"constant"` for population-invariant movement, or a list of integer
  vectors partitioning populations into blocks that share parameters.

- Movement_ageblk_spec:

  `"constant"` for age-invariant movement, or a list of integer vectors
  defining age blocks.

- Movement_yearblk_spec:

  `"constant"` for time-invariant movement, or a list of integer vectors
  defining year blocks.

- Movement_sexblk_spec:

  `"constant"` for sex-invariant movement, or a list of integer vectors
  defining sex blocks.

- Movement_seasblk_spec:

  `"constant"` for season-invariant movement, or a list of integer
  vectors defining season blocks.

- use_fixed_movement:

  Integer flag. `1` = movement rates are externally supplied; all
  `move_pars` are mapped to `NA` and not estimated. `0` = movement is
  estimated.

## Value

The input `input_list` with three `$map` entries updated:

- `$map$move_pars`:

  Factor vector for unstructured Markov transition parameters. Under
  `move_type = 0` with estimated movement, cells within the same block
  receive the same integer index; cells outside a spatial model or with
  fixed movement are `NA`. Entirely `NA` under `move_type = 1`.

- `$map$log_move_diffusion_pars`:

  Factor vector for CTMC diffusion parameters. Active (sequential
  integers) under `move_type = 1`; entirely `NA` under `move_type = 0`.

- `$map$move_preference_pars`:

  Factor vector for CTMC preference (taxis) parameters. Active under
  `move_type = 1`; entirely `NA` under `move_type = 0`.
