# Season weights for a fishery index reported once a year

A catch per unit effort reported for the year is the year's catch over
the year's effort. With fishing mortality proportional to effort, that
is each season's index weighted by the season's fully selected fishing
mortality, \\\sum\_\tau F\_\tau I\_\tau / \sum\_\tau F\_\tau\\, rather
than the seasons added together. A year with no fishing weighs its
seasons equally.

## Usage

``` r
get_fish_idx_seas_wt(F_seas)
```

## Arguments

- F_seas:

  Fully selected fishing mortality of one fleet in one region and year,
  one value per season.

## Value

Weights over the seasons, summing to one.
