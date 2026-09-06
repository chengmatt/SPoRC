# Draw state-space numbers-at-age innovations

Draws \\\eta\\ for a whole replicate at once. Drawing year by year would
only work when the year dim is independent or Markov; a separable
autoregression or a three-dimensional field correlates the whole span,
so the array is built up front and applied as the year loop reaches each
boundary.

## Usage

``` r
draw_naa_innovations(sim_env)
```

## Arguments

- sim_env:

  Simulation environment with the settings from
  [`Setup_Sim_NAA_state`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md)
  and the dimensions.

## Value

Array `[pop, region, year, season, age, sex]` of innovations, zero
outside the active ages, years and seasons.

## Details

Correlation is imposed dim by dim on independent normals, applying each
dim's Cholesky factor in turn, which is the reverse of how the penalty
whitens them. The three-dimensional field is the exception: its cohort
term couples age and year, so those two dims are drawn together from the
sparse precision rather than separately.
