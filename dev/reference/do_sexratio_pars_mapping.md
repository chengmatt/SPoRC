# Map sex ratio parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `sexratio_pars`, the proportion
of recruits assigned to the first sex. The `sexratio_pars` array has
dimensions `[n_pop x n_regions x n_sexratio_blocks]`, where
`n_sexratio_blocks` is the maximum number of time blocks across all
population-region combinations as defined in `$data$sexratio_blocks`.

## Usage

``` r
do_sexratio_pars_mapping(input_list, sexratio_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_sexes`, `$data$n_pop`, `$data$n_regions`,
  `$data$sexratio_blocks`, and `$data$rec_region_prop_spec` to be set by
  upstream setup functions.

- sexratio_spec:

  Character string specifying the estimation structure for
  `sexratio_pars`. One of:

  `"est_all"`

  :   Separate sex ratio parameter per population x region x block. Not
      permitted when `rec_region_prop_spec = 1` and `n_pop > 1`.

  `"est_shared_r"`

  :   Separate sex ratio per population x block, shared across regions
      within each population. Block membership is checked per region to
      ensure only valid blocks are assigned.

  `"est_shared_pop_r"`

  :   Single sex ratio per block, shared across all populations and
      regions. Requires identical block structures across all
      population-region combinations.

  `"fix"`

  :   All `sexratio_pars` fixed at their starting values (mapped to
      `NA`). Required when `n_sexes = 1`.

## Value

The input `input_list` with `$map$sexratio_pars` set to a factor vector
of length `prod(dim(par$sexratio_pars))`. Active parameters receive
sequential integer indices; fixed or invalid cells are `NA`.

## Details

When `n_sexes = 1`, estimation is meaningless and `sexratio_spec` must
be `"fix"`. Under the no-dispersal constraint
(`rec_region_prop_spec = 1` with `n_pop > 1`), `"est_all"` is prohibited
because non-natal regions receive no recruitment and cannot support
independent sex ratio estimates. The `"est_shared_pop_r"` option
additionally requires that all population-region combinations share the
same block structure; a mismatch raises an error.
