# Map stray rate parameters

Builds the factor map for `stray_rate_pars`
`[n_pop x max_stray_blocks]`, the logit-scale stray rates, where
`max_stray_blocks` is the largest number of time blocks across
populations in `$data$stray_rate_blocks`. Every parameter is fixed when
`n_pop = 1` or `use_fixed_stray_rate = 1`, whatever `stray_rate_spec`
asks for. Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

## Usage

``` r
do_stray_rate_mapping(input_list, stray_rate_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`. Requires `$data$n_pop`,
  `$data$stray_rate_blocks` and `$data$use_fixed_stray_rate`.

- stray_rate_spec:

  Estimation structure. `"fix"` holds every parameter at its starting
  value, `"est_all"` estimates per population and block, and
  `"est_shared_pop"` gives one parameter per block shared across
  populations, which requires identical block structures and errors
  otherwise.

## Value

`input_list` with `$map$stray_rate_pars` set to a factor vector of
length `prod(dim(par$stray_rate_pars))`. Active parameters take
sequential integers, fixed ones are `NA`.
