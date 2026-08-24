# Mean size of the plus group

Adjust for fish older than the accumulator age: the plus group is a
mixture of ages from the accumulator age onward, their numbers decaying
at \\Z = 0.2\\ per year, and their size growing linearly from the
curve's value at the accumulator age to \\L\_\infty\\ over a second
lifetime, \$\$\bar L\_+ = \frac{\sum\_{a=0}^{n} e^{-0.2 a}\\\[L_n +
(a/n)(L\_\infty - L_n)\]}{\sum\_{a=0}^{n} e^{-0.2 a}}\$\$ with \\n\\ the
accumulator age.

## Usage

``` r
plus_group_size(L_acc, Linf, n_acc)
```

## Arguments

- L_acc:

  Mean length at the accumulator age from the curve.

- Linf:

  Asymptotic length.

- n_acc:

  The accumulator age.

## Value

The adjusted plus-group mean length.
