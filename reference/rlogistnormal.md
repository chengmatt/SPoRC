# Simulate logistic normal variables

Simulate logistic normal variables

## Usage

``` r
rlogistnormal(exp, pars, comp_like, n_sexes)
```

## Arguments

- exp:

  Expected values

- pars:

  Parameters for a logistic normal (iid == 1 parameter, AR1 == 2
  parameters, 2D, by age and sex == 3 parameters, 3D, by age, sex, and
  region == 4 parameters)

- comp_like:

  Likelihood structure (iid == 2, ar1 == 3, 2d == 4, 3d == 5)

- n_sexes:

  Number of sexes
