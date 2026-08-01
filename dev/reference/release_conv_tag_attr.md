# Distribute Tagged Fish Releases to Full Population Dimensions

When tag release data are not recorded at full population resolution,
i.e. when one or more of the population, age, or sex dimensions are
unattended in `tag_attr`, this function distributes the known tag totals
to full `[n_pop, n_ages, n_sexes]` resolution using apportionment
weights derived from the release platform (population abundance, fishery
catch-at-age, or survey index-at-age). If all three dimensions are
attended (`tag_attr = "p_a_s"`), `tagged_fish` is returned unchanged
with no computation performed.

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

  Numeric vector or array of released tagged fish for a single tag
  cohort. Unattended dimensions are expected to be collapsed to index 1.
  Reshaped internally to `[n_pop, n_ages, n_sexes]`.

- tag_attr:

  Character string specifying which population dimensions are attended
  in `tagged_fish`. Constructed from any combination of `"p"`
  (population), `"a"` (age), and `"s"` (sex), joined by underscores
  (e.g. `"p_a_s"`, `"a"`, `"p_a"`). See
  [`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
  for the full set of valid strings.

- tag_release_platform:

  Character vector of length 2. Element 1 is the release platform: one
  of `"population"`, `"fishery"`, or `"survey"`. Element 2 is the fleet
  index as a character string (coerced to integer internally), or `NA`
  when `platform = "population"`.

- srv_sel:

  Numeric array of survey selectivity
  `[n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets]`. Used as the
  age-sex apportionment weight when `platform = "survey"`.

- fish_sel:

  Numeric array of fishery selectivity
  `[n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets]`. Used as the
  age-sex apportionment weight when `platform = "fishery"`.

- NAA:

  Numeric array of numbers-at-age prior to movement
  `[n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes]`. Used directly as
  weights when `platform = "population"`, and multiplied by selectivity
  for fishery and survey platforms.

- ty:

  Integer. Model year index of the tag release cohort.

- tseas:

  Integer. Season index of the tag release cohort.

- tr:

  Integer. Region index of the tag release cohort.

- n_pop:

  Integer. Number of populations.

- n_ages:

  Integer. Number of age classes.

- n_sexes:

  Integer. Number of sexes (1 or 2).

## Details

Apportionment weights are constructed from numbers-at-age
(`platform = "population"`), numbers-at-age multiplied by fishery
selectivity (`platform = "fishery"`), or numbers-at-age multiplied by
survey selectivity (`platform = "survey"`), all evaluated at the release
region, year, and season. Weights are then normalised conditionally on
the attended dimensions: the denominator for cell `[p, a, s]` is the sum
of raw weights across all cells that share the same indices in the
attended dimensions. This ensures that the marginal totals of
`tagged_fish` are preserved exactly along every attended dimension. For
example, if only age is attended (`tag_attr = "a"`), age-specific totals
in `tagged_fish` are preserved while tags are distributed across
population and sex in proportion to the platform weights.
