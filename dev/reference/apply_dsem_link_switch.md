# Take the linked deviations out of their own penalty

Every deviation penalty reads `data$map_<parameter>` rather than the map
itself, so zeroing out the linked cells there takes those deviations out
of it. Call from `sync_dev_map_data`, which rebuilds it before
`MakeADFun`.

## Usage

``` r
apply_dsem_link_switch(data)
```

## Arguments

- data:

  Data list holding `dsem_link_par` and `dsem_link_cell`.

## Value

`data` with the linked cells blanked in each mirror.
