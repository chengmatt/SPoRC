# Map Beverton-Holt steepness parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `steepness_h`, the
Beverton-Holt steepness parameter. The `steepness_h` array has
dimensions `[n_pop x n_regions]`.

## Usage

``` r
do_h_mapping(input_list, h_spec, rec_dd)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$rec_model`, `$data$n_pop`, `$data$n_regions`, and
  `$data$rec_dd` to be set by upstream setup functions.

- h_spec:

  Character string specifying the sharing structure for `steepness_h`,
  or `NULL` to estimate steepness independently across all relevant
  dimensions. Options when non-`NULL`:

  `"est_shared_pop_r"`

  :   Single steepness value shared across all populations and regions.
      All elements of `steepness_h` share factor level `1`. Required
      when `rec_dd = "global"` and `n_regions > 1`.

  `"est_shared_r"`

  :   Separate steepness per population, shared across regions within
      each population. Produces `n_pop` estimated parameters. Also valid
      under global density dependence.

  `"fix"`

  :   All `steepness_h` parameters fixed at their starting values
      (mapped to `NA`).

  `NULL`

  :   Estimate steepness independently: by population when `n_pop > 1`
      (shared across regions within each population), or by region when
      `n_pop = 1`. Not permitted when `rec_dd = "global"` and
      `n_regions > 1`.

- rec_dd:

  Recruitment density-dependence structure inherited from
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).
  `"global"` restricts valid `h_spec` values to `"est_shared_r"`,
  `"est_shared_pop_r"`, or `"fix"`, and prohibits `NULL`.

## Value

The input `input_list` with `$map$steepness_h` set to a factor vector of
length `prod(dim(par$steepness_h))`. Active parameters receive
sequential integer indices; unused parameters (mean recruitment model or
fixed steepness) are `NA`.

## Details

When `rec_model = 0` (mean recruitment), steepness has no role in the
stock-recruit relationship and all elements are mapped to `NA`
regardless of `h_spec`. For Beverton-Holt recruitment (`rec_model = 1`),
mapping follows `h_spec` subject to the density-dependence constraint:
global density dependence requires steepness to be shared across regions
(`"est_shared_r"` or `"est_shared_pop_r"`), since a single pooled
spawner-recruit relationship cannot support region-specific steepness
values.
