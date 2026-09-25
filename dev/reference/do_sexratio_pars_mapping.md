# Map sex ratio parameters

Builds the factor map for `sexratio_pars`
`[n_pop x n_regions x n_sexratio_blocks]`, the proportion of recruits
assigned to the first sex, where `n_sexratio_blocks` is the largest
number of time blocks across population and region in
`$data$sexratio_blocks`. Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

## Usage

``` r
do_sexratio_pars_mapping(input_list, sexratio_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`. Requires `$data$n_sexes`,
  `$data$n_pop`, `$data$n_regions`, `$data$sexratio_blocks` and
  `$data$rec_region_prop_spec`.

- sexratio_spec:

  Estimation structure for `sexratio_pars`. `"est_all"` gives one
  parameter per population, region and block, and is not permitted when
  `rec_region_prop_spec = 1` and `n_pop > 1`. `"est_shared_r"` gives one
  per population and block, shared across its regions.
  `"est_shared_pop_r"` gives one per block across every population and
  region, and requires identical block structures. `"fix"` holds every
  parameter at its starting value and is required when `n_sexes = 1`.

## Value

`input_list` with `$map$sexratio_pars` set to a factor vector of length
`prod(dim(par$sexratio_pars))`. Active parameters take sequential
integers, fixed or invalid cells are `NA`.
