# Preference terms that vary over years

Internal helper called by
[`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md).
A dsem on the movement deviations describes the year to year part of
habitat preference, so a preference covariate that also varies over
years writes that part a second time, and the formula's version of it is
not penalized.

## Usage

``` r
get_yr_varying_pref_terms(input_list)
```

## Arguments

- input_list:

  Named list with `$data`, built through
  [`Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md).

## Value

Character vector naming the preference design columns whose values
change across years within a population, region, season, age and sex.
Empty when movement is not CTMC, when the preference formula has no
terms, or when none of its terms vary over years.
