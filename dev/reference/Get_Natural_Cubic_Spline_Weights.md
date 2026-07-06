# Construct a natural cubic spline interpolation weight matrix

Precomputes a linear operator \\W\\ such that `W %*% y_nodes` reproduces
natural cubic spline interpolation (zero second derivative at the
endpoints) of the node values `y_nodes` (defined at `x_nodes`) onto the
query points `x_out`. Because natural cubic spline interpolation is
linear in the node values for fixed node/query positions, \\W\\ depends
only on `x_nodes` and `x_out` (both treated as fixed data), never on the
node values themselves. This lets bicubic/cubic selectivity splines (see
[`Get_Selex`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex.md),
`Selex_Model == 8`) be evaluated as a pair of matrix multiplications
against AD parameter vectors, rather than re-solving a spline system on
every function evaluation.

## Usage

``` r
Get_Natural_Cubic_Spline_Weights(x_nodes, x_out)
```

## Arguments

- x_nodes:

  Numeric vector of strictly increasing node (knot) positions, length
  \\n \ge 1\\. Typically bin or year positions rescaled to \\\[0,1\]\\.

- x_out:

  Numeric vector of query positions at which the spline is to be
  evaluated, length \\m\\. Values are clamped to the innermost spline
  segment if they fall outside `range(x_nodes)` (no extrapolation).

## Value

Numeric \\m \times n\\ weight matrix \\W\\. When \\n == 1\\ (a single
node), \\W\\ is a column of ones (the interpolated curve is constant,
equal to the single node value). When \\n == 2\\, natural boundary
conditions force the spline to reduce to linear interpolation.

## Details

Natural cubic spline second derivatives \\M\\ at the nodes solve a
tridiagonal linear system \\A M = R y\\ where \\A\\ and \\R\\ depend
only on the node spacing `diff(x_nodes)` (data), so \\M = A^{-1} R y =
\text{Mmat} \\ y\\ is itself linear in \\y\\. Substituting into the
standard piecewise-cubic evaluation formula for each query point yields
one row of \\W\\ per query point, each a fixed linear combination of the
node basis vectors and rows of `Mmat`.

## Examples

``` r
if (FALSE) { # \dontrun{
x_nodes <- seq(0, 1, length.out = 4)
x_out <- seq(0, 1, length.out = 20)
W <- Get_Natural_Cubic_Spline_Weights(x_nodes, x_out)
y_nodes <- c(0.1, 0.8, 0.6, 0.3)
W %*% y_nodes # interpolated curve at x_out
} # }
```
