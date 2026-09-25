# Draw covariate observations from their family and link

One draw per cell: the state through the link is the mean, and the
family draws about it, as
[`get_dsem_obs_nLL`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_obs_nLL.md)
evaluates it. Fixed returns the state itself.

## Usage

``` r
draw_dsem_cov_obs(state, family, link, obs_sd, tweedie_p, fixed_sd = NULL)
```

## Arguments

- state:

  Array `[year, 1, sim]` of link-scale cells.

- family:

  Family code as in
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md).

- link:

  Link code.

- obs_sd:

  Spread parameter on the family's own scale.

- tweedie_p:

  Tweedie power in (1, 2).

- fixed_sd:

  Known sd per year for gaussian_fixed_sd, recycled over replicates;
  NULL otherwise.

## Value

Array shaped like `state`.
