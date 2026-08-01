# Map stray rate parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `stray_rate_pars`, the
logit-scale stray rate parameters. The array has dimensions
`[n_pop x max_stray_blocks]`, where `max_stray_blocks` is the maximum
number of time blocks across all populations as defined by
`$data$stray_rate_blocks`.

## Usage

``` r
do_stray_rate_mapping(input_list, stray_rate_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_pop`, `$data$stray_rate_blocks`, and
  `$data$use_fixed_stray_rate` to be set by upstream functions.

- stray_rate_spec:

  Character string specifying the estimation structure. One of:

  `"fix"`

  :   All parameters fixed at starting values (mapped to `NA`).

  `"est_all"`

  :   Estimate independently per population x block. Produces
      `n_pop x n_unique_blocks` estimated parameters.

  `"est_shared_pop"`

  :   Single parameter per block, shared across all populations.
      Requires identical block structures across all populations. An
      error is raised if block indices differ.

## Value

The input `input_list` with `$map$stray_rate_pars` set to a factor
vector of length `prod(dim(par$stray_rate_pars))`. Active parameters
receive sequential integer indices; fixed parameters are `NA`.

## Details

When `n_pop = 1`, straying is not applicable and all parameters are
automatically fixed to `NA`. When `use_fixed_stray_rate = 1`, the
objective function reads from `fixed_stray_rate` directly and
`stray_rate_pars` are not used; all elements are fixed regardless of
`stray_rate_spec`.
