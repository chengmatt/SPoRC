# Run Iterative Francis Reweighting Procedure

Runs an iterative Francis reweighting procedure for composition data
(fishery and survey age- and length-compositions). The function
reweights input data, repeatedly fits the model, and computes updated
Francis weights.

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

  A list of model input data, including at least observed compositions
  (\`ObsFishAgeComps\`, \`ObsFishLenComps\`, \`ObsSrvAgeComps\`,
  \`ObsSrvLenComps\`) and corresponding weights (\`Wt_FishAgeComps\`,
  \`Wt_FishLenComps\`, \`Wt_SrvAgeComps\`, \`Wt_SrvLenComps\`).

- parameters:

  A list of model parameters to be passed to \[fit_model()\].

- mapping:

  A list or mapping object used to specify fixed or estimated parameters
  in \[fit_model()\].

- random:

  A character string of random effects passed to \[fit_model()\].

- n_francis_iter:

  Integer. Number of Francis reweighting iterations to perform. Default
  is \`10\`.

- newton_loops:

  Integer. Number of Newton loops passed to \[fit_model()\]. Default is
  \`0\`.

## Value

A list with three elements:

- obj:

  The fitted model object returned by \[fit_model()\], including all
  elements of a TMB object, data, parameters, mapping, random effects
  specified, and report.

- end_mean_francis:

  A summary of the mean Francis weights from the first iteration.

- end_mean_francis:

  A summary of the mean Francis weights from the final iteration.

- recorded_weights:

  A summary of recorded francis weights from each iteartion.

## See also

Other Francis Reweighting:
[`do_francis_reweighting()`](https://chengmatt.github.io/SPoRC/reference/do_francis_reweighting.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  out <- run_francis(data = data,
                     parameters = parameters,
                     mapping = mapping,
                     random = NULL,
                     n_francis_iter = 5,
                     newton_loops = 3)
  out$obj
  out$mean_francis
} # }
```
