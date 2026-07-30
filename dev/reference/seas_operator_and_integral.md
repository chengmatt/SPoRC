# Seasonal transition operator and abundance integral from one matrix exponential

Returns both the seasonal transition operator and the season-integrated
abundance operator for `move_timing = 2`, at the cost of a single matrix
exponential.

## Usage

``` r
seas_operator_and_integral(Z, Q, dur = 1)
```

## Arguments

- Z:

  Numeric vector of length `n_regions` giving total mortality for this
  season, already scaled by season duration (i.e. the same quantity
  stored in `ZAA`). Pass zeros for a movement-only step.

- Q:

  Square `[n_regions x n_regions]` instantaneous rate matrix (generator)
  in row convention, as stored in `Mrate`. Required when
  `move_timing = 2`; ignored otherwise.

- dur:

  Season duration used to scale `Q`. Should be `seasdur[seas]`. Only
  used when `move_timing = 2`.

## Value

A list with

- `T`:

  The seasonal transition operator in row convention, identical to
  `build_seas_operator(..., move_timing = 2)`.

- `Integral`:

  \\\int_0^1 e^{A\tau}d\tau\\ in column convention; `Integral %*% N` is
  the season-integrated abundance.

## Details

The Van Loan block used to recover the integral already contains the
transition operator in its top-left corner,
\$\$\exp\left(\begin{bmatrix} A & I \\ 0 & 0\end{bmatrix}\right) =
\begin{bmatrix} e^{A} & \int_0^1 e^{A\tau}d\tau \\ 0 &
I\end{bmatrix},\$\$ so callers that need the population step *and* the
catch integral for the same stratum – which is every fished stratum,
since the dynamics advance the numbers and the Baranov equation
integrates them over the identical \\A\\ – should take both from here
rather than exponentiating \\A\\ once and the block again. Under
reverse-mode AD the adjoint of a matrix exponential is far more
expensive than its forward evaluation, so halving the number of
exponentials on the tape is worth more than the flop count alone
suggests.
