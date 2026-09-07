# Map AR1 correlation parameter for recruitment deviations

Constructs the `RecDevs_rho` factor map. `RecDevs_rho` is only read when
`RecDevs_model = "ar1"` (see
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md));
under any other `RecDevs_model` every `RecDevs_rho` parameter is mapped
to `NA` regardless of `RecDevs_rho_spec`, since the recruitment penalty
never reads it.

## Usage

``` r
do_RecDevs_rho_mapping(input_list, RecDevs_rho_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions.

- RecDevs_rho_spec:

  Character string controlling the sharing and estimation structure for
  `RecDevs_rho`: one of `"est_all"`, `"est_shared_pop"`,
  `"est_shared_r"`, `"est_shared_pop_r"`, or `"fix"`.

## Value

The input `input_list` with `$map$RecDevs_rho` set to a factor vector of
length `prod(dim(par$RecDevs_rho))`.
