# Classify conventional tag likelihood family

Maps a conventional tag likelihood type code to a coarse family label
used for OSA packing and evaluation.

## Usage

``` r
tag_fam_of(lt)
```

## Arguments

- lt:

  Integer likelihood type code:

  - 0, 1: count-based (Poisson / NB)

  - 2, 3, 4, 5: composition-based (multinomial / Dirichlet-multinomial)

## Value

A character scalar: \`"count"\`, \`"comp"\`, or \`NA_character\_\` if
the code is not recognized.
