# Put an OSA observation slice into the order it is peeled in

[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
predicts observations in the order its `ord` attribute gives rather than
in storage order, and the conditional factorizations below are the right
ones only when they follow that order. RTMB's own `dmultinom` OSA method
and WHAM's `age_comp_osa.hpp` apply the same permutation. Under the
ascending subset
[`osa_keep_subset`](https://chengmatt.github.io/SPoRC/dev/reference/osa_keep_subset.md)
builds the two orders already agree, so this is the identity on every
call the package makes; it is what keeps a caller-supplied `subset`
correct.

## Usage

``` r
osa_order(xobs, n)
```

## Arguments

- xobs:

  An `"osa"` object or plain numeric vector.

- n:

  Length of the slice.

## Value

Integer permutation, the identity when there is nothing to reorder.
