# Integrate abundance over a season under continuous movement

Computes \\\int_0^1 e^{A\tau} d\tau \\ N(0)\\, the season-integrated
abundance required for the spatial Baranov catch equation when movement
and mortality act simultaneously. Under `move_timing = 2` the usual
\\F/Z \\ (1 - e^{-Z}) N\\ form is not valid, because fish redistribute
among regions while they are dying.

## Usage

``` r
integrate_seas_abundance(N, Z, Q, dur = 1)
```

## Arguments

- N:

  Numeric vector of length `n_regions` at the start of the season.

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

Numeric vector of length `n_regions` giving season-integrated abundance
by region. Multiplying elementwise by a fishing mortality vector gives
catch.

## Details

The integral runs over the unit interval, not over \\\[0, \Delta\]\\,
because \\A = Q^\top \Delta - \mathrm{diag}(Z)\\ already carries the
season duration in both of its terms. The integration variable is
elapsed *fraction* of the season, so a full season is \\\tau = 1\\.

Van Loan's block identity \$\$\exp\left(\begin{bmatrix} A & I \\ 0 &
0\end{bmatrix}\right) = \begin{bmatrix} e^{A} & \int_0^1 e^{A\tau}d\tau
\\ 0 & I\end{bmatrix}\$\$ which recovers the integral from a single
matrix exponential of twice the dimension. This is preferred over
\\A^{-1}(e^{A} - I)\\ because it avoids an explicit inverse and stays
well conditioned as \\A\\ approaches singularity. That is not an edge
case here: the generator has a zero eigenvalue (its rows sum to zero),
so \\A\\ becomes exactly singular as total mortality approaches zero,
precisely where the integral itself remains perfectly well behaved.
