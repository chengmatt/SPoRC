# Get Design Matrices for CTMC Movement

Constructs the design matrices for the diffusion and preference
components of a Continuous Time Markov Chain (CTMC) movement model.
These matrices are used to parameterize movement rates in terms of
covariates specified by formulas.

## Usage

``` r
get_movement_dp_design_matrix(data, preference_formula, diffusion_formula)
```

## Arguments

- data:

  A `data.frame` containing the CTMC covariates. Must include all
  variables referenced in `diffusion_formula` and `preference_formula`.

- preference_formula:

  An R formula describing the linear predictor for movement preference
  (taxis). Variables must exist in `data`.

- diffusion_formula:

  An R formula describing the linear predictor for diffusion rates.
  Variables must exist in `data`.

## Value

A `list` with the following components:

- `n_theta`:

  Number of diffusion parameters (columns in `X_zk`).

- `n_gamma`:

  Number of preference parameters (columns in `W_zk`).

- `X_zk`:

  Diffusion design matrix constructed from `diffusion_formula` and
  `data`.

- `W_zk`:

  Preference (taxis) design matrix constructed from `preference_formula`
  and `data`.
