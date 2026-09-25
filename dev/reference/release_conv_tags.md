# Release conventional tags in the simulation

Spreads each cohort's `n_tags`, or `n_tags_rel_input` when supplied,
across populations, ages and sexes in proportion to the
selectivity-weighted abundance `NAA_bef` of the release platform,
rounding to integers. `conv_fish_tag_attr` is then applied through
[`marginalize_conv_fish_tags`](https://chengmatt.github.io/SPoRC/dev/reference/marginalize_conv_fish_tags.md)
to give the observation-level release array `conv_tagged_fish_attr`, at
the resolution the recapture likelihood reads.

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

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place: `$conv_tagged_fish` and `$conv_tagged_fish_attr`
  for each cohort released in year `y`.

## Value

`invisible(NULL)`; everything is modified by reference within `sim_env`.

## Details

Without `n_tags_rel_input`, a survey or fishery platform scales the
region's total against the selectivity-weighted global abundance, and a
population platform against the region's share of total abundance.
