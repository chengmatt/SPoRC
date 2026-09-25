# Map Beverton-Holt steepness parameters

Builds the factor map for `steepness_h` `[n_pop x n_regions]`. Every
element is `NA` under `rec_model = 0`, where steepness has no role.
Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

## Usage

``` r
do_h_mapping(input_list, h_spec, rec_dd)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`. Requires
  `$data$rec_model`, `$data$n_pop`, `$data$n_regions` and
  `$data$rec_dd`.

- h_spec:

  Sharing structure for `steepness_h`. `"est_shared_pop_r"` gives one
  value across every population and region, required when
  `rec_dd = "global"` and `n_regions > 1`. `"est_shared_r"` gives one
  per population. `"fix"` holds every value at its starting value.
  `NULL` estimates by population when `n_pop > 1` and by region when
  `n_pop = 1`, and is not permitted under global density dependence with
  `n_regions > 1`.

- rec_dd:

  Density dependence inherited from
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).
  `"global"` restricts `h_spec` to `"est_shared_r"`,
  `"est_shared_pop_r"` or `"fix"`.

## Value

`input_list` with `$map$steepness_h` set to a factor vector of length
`prod(dim(par$steepness_h))`. Active parameters take sequential
integers, unused ones are `NA`.
