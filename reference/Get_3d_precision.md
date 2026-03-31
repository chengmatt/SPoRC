# Constructor algorithin for correlations within ages, years, and cohort (Work in Progress)

Constructor algorithin for correlations within ages, years, and cohort
(Work in Progress)

## Usage

``` r
Get_3d_precision(
  n_ages,
  n_yrs,
  pcorr_age,
  pcorr_year,
  pcorr_cohort,
  ln_var_value,
  Var_Type
)
```

## Arguments

- n_ages:

  Number of ages

- n_yrs:

  Number of years

- pcorr_age:

  correlations for age

- pcorr_year:

  correaltions for year

- pcorr_cohort:

  correlaitons for cohort

- ln_var_value:

  log space variance

- Var_Type:

  variance type == 0, marginal (stationary and slower run time), == 1
  conditional (non-statationary, faster run time)

## Value

Sparse precision matrix dimensioned by n_ages \* n_years, n_ages \*
n_years
