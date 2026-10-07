# Warn when an aggregated composition is flagged outside region one

An aggregated composition (`"agg"`) is one over every region and sex,
kept in region one, where the default input sample size is written and
where the operating model draws it. The likelihood reads the first
region flagged, so a flag in another region is a region-resolved
composition fit against the whole model's, or a second copy of the same
one that is never read.

## Usage

``` r
check_agg_comp_regions(use_arr, comp_type, arg_name)
```

## Arguments

- use_arr:

  Use array, region by year by season by fleet, with a leading
  population dim for a population-specific data source.

- comp_type:

  Composition type matrix `[n_years, n_fleets]`, 0 for aggregated.

- arg_name:

  Name of the `Use` argument, used in the warning.

## Value

`NULL`, invisibly. Called for the warning it raises.
