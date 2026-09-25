# Series a dsem can link, and where each one's cells sit

One series per deviation array, kept at every index except the year:
`ln_growth_devs` of a two sex model gives one series per population,
region, parameter and sex. An arrow names the ones to link.

## Usage

``` r
get_dsem_link(input_list, dsem_processes, arrow_text, n_grid_yrs)
```

## Arguments

- input_list:

  List with `data`, `par` and `map`.

- dsem_processes:

  Labels from `dsem_process_table` to offer.

- arrow_text:

  The arrow lines as one string, comments already removed.

- n_grid_yrs:

  Rows of the dsem grid.

## Value

List with `offered` (every name the processes could give, with
`offered_label`, `offered_par` and `offered_cell`) and, for the linked
ones, `name`, `par`, `label`, `sigma_par`, `penalty_reads_map`,
`grid_row` (the grid rows that array reaches) and `cell` (where each of
those rows sits in the flattened array).
