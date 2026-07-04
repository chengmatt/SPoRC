# Marginalise conventional fishery tag arrays across unattended dimensions

Collapses population, age, and/or sex dimensions of a tag count array
`[n_pop × n_ages × n_sexes]` by summing over dimensions absent from
`tag_recaptures_attr`, placing the result into index 1 of the
corresponding dimension and zeroing all other indices. Region and fleet
are not handled here (they are managed at the calling level). The
function is used to align the release array `conv_tagged_fish` with the
attended resolution of the recapture likelihood.

## Usage

``` r
marginalize_conv_fish_tags(vals, tag_recaptures_attr, n_pop, n_ages, n_sexes)
```

## Arguments

- vals:

  Numeric vector or array of tag counts, interpreted as a
  `[n_pop × n_ages × n_sexes]` array.

- tag_recaptures_attr:

  Character string specifying attended dimensions. Same format as
  `conv_fish_tag_attr` in
  [`Setup_Sim_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md):
  any combination of `"p"`, `"a"`, `"s"` joined by underscores.

- n_pop:

  Integer. Number of populations.

- n_ages:

  Integer. Number of age classes.

- n_sexes:

  Integer. Number of sexes.

## Value

Array `[n_pop × n_ages × n_sexes]` with unattended dimensions summed
into index 1 and all other indices set to zero.
