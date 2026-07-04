# Release conventional tags in the simulation

For each tag cohort scheduled for release in year `y`, distributes
`n_tags` (or `n_tags_rel_input` if provided) across populations, ages,
and sexes proportional to the selectivity-weighted abundance (`NAA_bef`)
of the release platform (survey, fishery, or population). Tagged fish
counts are rounded to integers. The attended attribute string
`conv_fish_tag_attr` is then applied via
[`marginalize_conv_fish_tags`](https://chengmatt.github.io/SPoRC/dev/reference/marginalize_conv_fish_tags.md)
to produce the observation-level release array `conv_tagged_fish_attr`,
which is consistent with the dimension resolution of the recapture
likelihood.

## Usage

``` r
release_conv_tags(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment created by
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).
  Modified in place: `$conv_tagged_fish[tc, , , , sim]` and
  `$conv_tagged_fish_attr[tc, , , , sim]` for each cohort released in
  year `y`.

## Value

`invisible(NULL)`. All modifications are made by reference within
`sim_env`.

## Details

For survey and fishery platforms, total tags in the release region are
scaled relative to the selectivity-weighted global abundance to allocate
region-specific cohort sizes when `n_tags_rel_input` is not provided.
For the population platform, scaling is proportional to the region's
share of total abundance.
