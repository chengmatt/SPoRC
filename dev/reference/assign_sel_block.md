# Assign a value to every region x year cell of one fleet belonging to a selectivity block

Selectivity block arrays are `[region, year, fleet]`, so one fleet's
slice is a region by year matrix. `which(slice == block)` on it counts
down the columns, giving positions in `1:(n_regions * n_years)` rather
than years. Used as a year subscript those are wrong whenever there is
more than one region: either the subscript is out of bounds, or, when
the block is early enough that the positions stay below `n_years`, the
wrong years are written with no error at all. With three regions and 35
years, a block covering years 1-5 writes years 1-15. With one region the
position equals the year, which is why this only shows up in spatial
models.

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
