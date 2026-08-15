# Assign a value to every region x year cell of one fleet belonging to a selectivity block

Selectivity block arrays are `[region, year, fleet]`, so a single
fleet's slice is a `region x year` MATRIX. `which(slice == block)` on
that matrix returns LINEAR positions running down the columns, in
`1:(n_regions * n_years)` – not year indices. Using them as a year
subscript (`arr[, which(...), fleet] <- value`) is therefore wrong
whenever `n_regions > 1`: it either errors with a subscript out of
bounds, or, when the block is early enough that the linear positions
stay within `n_years`, SILENTLY writes the wrong years. With three
regions and 35 years, a block covering years 1-5 produces linear
positions 1-15 and quietly overwrites years 1-15. At `n_regions == 1`
the linear position equals the column index, which is why this only
shows up in spatial models.

## Usage

``` r
assign_sel_block(arr, blocks_arr, fleet, block, value)
```

## Arguments

- arr:

  Array `[region, year, fleet]` to write into.

- blocks_arr:

  Block array `[region, year, fleet]`, same first three dims as `arr`.

- fleet:

  Fleet index.

- block:

  Block value to match.

- value:

  Scalar to assign to the matching cells.

## Value

`arr` with the matching cells of `fleet` set to `value`.

## Details

Indexing with the logical matrix directly is correct in both cases, and
stays correct if blocks are ever allowed to differ between regions.
