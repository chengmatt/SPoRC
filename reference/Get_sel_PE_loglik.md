# Get selectivity process error log likelihoods (positive)

Get selectivity process error log likelihoods (positive)

## Usage

``` r
Get_sel_PE_loglik(
  PE_model,
  PE_pars,
  ln_devs,
  map_sel_devs,
  sel_vals,
  do_sel_pen
)
```

## Arguments

- PE_model:

  Process error model values (1, 2, 3, 4, and 5) (iid, random walk, 3d
  marginal, 3d conditional, and 2dar1)

- PE_pars:

  Process error parameters

- ln_devs:

  Deviations

- map_sel_devs:

  selectivity deviations to share

- sel_vals:

  Selectivity values (either length or age based)

- do_sel_pen:

  Whether temporal or bin penalties are used

## Value

numeric value of log likelihood (in positive space)
