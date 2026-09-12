# Refuse an at-age observation that no model age is read as

An observed age whose column of the fleet's ageing error is all zero is
predicted as zero whatever the population does, so an observation there
cannot be fit. This names the first fleet, year and observed age where
that happens.

## Usage

``` r
check_at_age_read(input_list, use, fleet_field, what, pop = FALSE)
```

## Arguments

- input_list:

  Named list with `$data`, after
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

- use:

  Use array of the data source.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- what:

  Name used in the message.

- pop:

  Logical. `TRUE` for the population-specific data source.

## Value

`invisible(NULL)`. Called for its error.
