# Run Iterative Francis Reweighting Procedure

Runs an iterative Francis reweighting procedure for pooled and
population-specific fishery and survey age and length compositions.
Repeatedly fits the model and updates composition weights until the
specified number of iterations is reached.

## Usage

``` r
run_francis(
  data,
  parameters,
  mapping,
  random = NULL,
  n_francis_iter = 10,
  newton_loops = 0
)
```

## Arguments

- data:

  A list of model input data containing observed compositions, usage
  flags, input sample sizes, and weight arrays for both pooled
  (`Wt_FishAgeComps`, `Wt_FishLenComps`, `Wt_SrvAgeComps`,
  `Wt_SrvLenComps`, `Wt_FishAgeComps_discard`,
  `Wt_FishLenComps_discard`) and population-specific
  (`Wt_FishAgeComps_pop`, `Wt_FishLenComps_pop`, `Wt_SrvAgeComps_pop`,
  `Wt_SrvLenComps_pop`, `Wt_FishAgeComps_discard_pop`,
  `Wt_FishLenComps_discard_pop`) data streams.

- parameters:

  A list of model parameters passed to
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- mapping:

  A list or mapping object passed to
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- random:

  Character vector of random effects passed to
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).
  Default `NULL`.

- n_francis_iter:

  Integer. Number of Francis reweighting iterations. Default `10`.

- newton_loops:

  Integer. Number of Newton refinement steps passed to
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).
  Default `0`.

## Value

A named list with:

- `obj`:

  The fitted model object from the final iteration, augmented with
  `$data`, `$parameters`, `$mapping`, `$random`, and `$rep`.

- `start_mean_francis`:

  Mean composition fits from the first iteration, as returned by
  [`do_francis_reweighting`](https://chengmatt.github.io/SPoRC/dev/reference/do_francis_reweighting.md).

- `end_mean_francis`:

  Mean composition fits from the final iteration.

- `recorded_weights`:

  Long-format dataframe of Francis weights from every iteration for all
  pooled and population-specific data streams, with columns `Region`,
  `Year`, `Seas`, `Sex`, `Fleet`, `Weight`, `Type`, `Pop` (pooled rows
  have `Pop = NA`), and `iter`.

## See also

Other Francis Reweighting:
[`do_francis_reweighting()`](https://chengmatt.github.io/SPoRC/dev/reference/do_francis_reweighting.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  out <- run_francis(data = data, parameters = parameters,
                     mapping = mapping, random = NULL,
                     n_francis_iter = 5, newton_loops = 3)
  out$obj
  out$end_mean_francis
  out$recorded_weights
} # }
```
