# Map AR1 correlation parameter for fishing mortality deviations

Constructs the `Fdev_rho` factor map. `Fdev_rho` is only meaningful when
`Fdev_model = "ar1"` (see
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md));
for any other `Fdev_model`, all `Fdev_rho` parameters are mapped to `NA`
regardless of `Fdev_rho_spec`, since they are unused by
[`Get_Fdev_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Fdev_PE_loglik.md).

## Usage

``` r
do_Fdev_rho_mapping(input_list, Fdev_rho_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions.

- Fdev_rho_spec:

  Character string controlling the sharing and estimation structure for
  `Fdev_rho`, following the same convention as
  [`do_sigmaF_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sigmaF_mapping.md)'s
  `sigmaF_spec`: one of `"est_all"`, `"est_shared_r"`,
  `"est_shared_seas"`, `"est_shared_f"`, `"est_shared_r_seas"`,
  `"est_shared_r_f"`, `"est_shared_seas_f"`, `"est_shared_r_seas_f"`, or
  `"fix"`.

## Value

The input `input_list` with `$map$Fdev_rho` set to a factor vector of
length `prod(dim(par$Fdev_rho))`.
