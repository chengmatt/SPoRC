# Calculate the Corrected marginal AIC (AICc) from Optimization Results

Computes the corrected marginal Akaike Information Criterion (AICc) for
model selection using optimization results. It supports objects returned
from different optimizers, such as \`optim\` or \`nlminb\`.

## Usage

``` r
marg_AIC(opt, p = 2, n = Inf)
```

## Arguments

- opt:

  A list containing optimization results. Must include either:

  - \`"par"\` and \`"objective"\` (e.g., from \`optim\`), or

  - \`"par"\` and \`"value"\` (e.g., from \`nlminb\`)

- p:

  Numeric. Penalty multiplier for the number of parameters. Default is
  2.

- n:

  Numeric. Sample size. Default is \`Inf\`.

## Value

Numeric. The corrected AIC (AICc) value.
