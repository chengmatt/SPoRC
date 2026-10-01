# Blocks of a movement deviation dim

Turns one switch of
[`Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md)
into a block id per level: `NA` for an inactive level, one block for
`"none"`, a block per active level for `"iid"`, `"ar1"` and `"us"`, and
the blocks as given for a list. A dsem splits every active level, since
it links one series per cell.

## Usage

``` r
move_dim_blocks(switch, active, n_levels, dsem = FALSE)
```

## Arguments

- switch:

  The switch's value.

- active:

  Integer vector of the active levels.

- n_levels:

  Number of levels of the dim.

- dsem:

  Logical, whether a dsem holds the deviations.

## Value

Integer vector of length `n_levels`, the block of each level.
