# Distribute Tagged Fish Releases to Full Population Dimensions

Spreads tag totals recorded at less than full resolution out to
`[n_pop, n_ages, n_sexes]`, using weights from the release platform.
`tagged_fish` comes back unchanged when every dim is attended
(`tag_attr = "p_a_s"`).

## Usage

``` r
release_conv_tag_attr(
  tagged_fish,
  tag_attr,
  tag_release_platform,
  srv_sel,
  fish_sel,
  NAA,
  ty,
  tseas,
  tr,
  n_pop,
  n_ages,
  n_sexes
)
```

## Arguments

- tagged_fish:

  Released tagged fish for one cohort, with the unattended dims
  collapsed to index 1 and reshaped internally to
  `[n_pop, n_ages, n_sexes]`.

- tag_attr:

  Which dims are attended, built from `"p"`, `"a"` and `"s"` joined by
  underscores. See
  [`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md).

- tag_release_platform:

  Character vector of length 2: the platform, `"population"`,
  `"fishery"` or `"survey"`, and the fleet index as a string, or `NA`
  under `"population"`.

- srv_sel, fish_sel:

  Survey and fishery selectivity
  `[n_regions, n_yrs, n_ages, n_sexes, n_fleets]`, used as the age and
  sex weight under their own platform.

- NAA:

  Numbers at age before movement
  `[n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes]`, the weights
  themselves under `"population"` and multiplied by selectivity
  otherwise.

- ty, tseas, tr:

  Year, season and region indices of the release cohort.

- n_pop, n_ages, n_sexes:

  Dimension sizes.

## Details

The weights are the numbers at age under `"population"`, or the numbers
at age times fishery or survey selectivity under the other two
platforms, all read at the release region, year and season. They are
normalized conditionally on the attended dims: the denominator for cell
`[p, a, s]` sums the raw weights over every cell sharing its indices in
those dims, so the marginal totals of `tagged_fish` are preserved
exactly along each of them. With only age attended, for instance, the
age totals are kept and the tags are spread over population and sex in
proportion to the weights.
