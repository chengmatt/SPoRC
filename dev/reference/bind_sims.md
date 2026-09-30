# Stack one array per replicate

The operating model reads every input with the replicates on the last
dim.

## Usage

``` r
bind_sims(parts)
```

## Arguments

- parts:

  List of arrays, one per replicate, all the same shape.

## Value

One array, replicates on the last dim.
