# Check every parameter block against the map that indexes it

RTMB pairs a parameter with its map by position, so the two must be the
same length. They come apart when a starting value is supplied at the
wrong shape, or when a map is built from dimensions the parameter does
not have. Neither is caught where it happens: the objective reads the
shorter of the two past its end, and RTMB reports an invalid advector
from somewhere unrelated.

## Usage

``` r
check_par_map_lengths(parameters, mapping)
```

## Arguments

- parameters:

  Parameter list.

- mapping:

  Map list. Entries naming a parameter that is absent are reported too,
  since a map with no parameter is silently ignored.

## Value

`invisible(NULL)`. Called for its error.
