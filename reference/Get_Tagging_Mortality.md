# Get values for total tag mortality and tag fishing mortality

Get values for total tag mortality and tag fishing mortality

## Usage

``` r
Get_Tagging_Mortality(
  tag_selex,
  tag_natmort,
  Fmort,
  natmort,
  Tag_Shed,
  fish_sel,
  n_regions,
  n_ages,
  n_sexes,
  n_fish_fleets,
  y,
  what
)
```

## Arguments

- tag_selex:

  Tag selectivity options, == 0 (uniform, with F from fleet 1 (dominant
  fleet)), == 1 (sex-averaged selectivity, with F from fleet 1 (dominant
  fleet)), == 2 (sex-specific selectivity, with F from fleet 1 (dominant
  fleet)), 3 (uniform with F summed across fleets), 4 (sex averaged with
  F summed across fleets, weighted sum), 5 (sex-specific with F summed
  across fleets, weighted sum)

- tag_natmort:

  Tag natural mortality options == 0 (averaged across sexes and ages),
  == 1 (averaged across sexes, but unique for ages), == 2 (sex-specific,
  but averaged across ages), == 3 (sex-and age-specific)

- Fmort:

  Array of fishing mortality, dimensioned by n_region, n_years,
  n_fish_fleets

- natmort:

  Array of fishing mortality, dimensioned by n_region, n_years, n_ages,
  n_sexes

- Tag_Shed:

  Scalar chronic tag shedding rate

- fish_sel:

  Array of fishery selectivity, dimensioned by n_region, n_years,
  n_ages, n_sexes, n_fish_fleetss

- n_regions:

  Number of regions

- n_ages:

  Number of ages

- n_sexes:

  Number of sexes

- n_fish_fleets:

  Number of fishery fleets

- y:

  Year index

- what:

  Whether to return Z or F (total or fishing mortality)

## Value

Z or F values from tagging specifications
