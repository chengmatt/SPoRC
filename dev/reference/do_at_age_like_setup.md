# Set the error scale and likelihood of one at-age stream

An at-age observation may be lognormal or normal, and its standard
deviation may come from an estimated parameter, from reported standard
errors, or from both. This is the parity the aggregated index streams
already have, stated per fleet.

## Usage

``` r
do_at_age_like_setup(
  input_list,
  like_type,
  sigma_form,
  stream,
  fleet_field,
  pop = FALSE
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- like_type:

  `"lognormal"` or `"normal"`, one setting for every fleet or one per
  fleet.

- sigma_form:

  `"none"` for the parameter alone, `"data"` for the reported standard
  errors alone, `"est_additive"` or `"est_quadrature"` for both.

- stream, fleet_field, pop:

  See
  [`do_at_age_type_setup`](https://chengmatt.github.io/SPoRC/dev/reference/do_at_age_type_setup.md).

## Value

`input_list` with `$data$<stream>_LikeType` and
`$data$<stream>_sigma_form` set.
