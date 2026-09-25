# Compute Movement Process Error Log-Likelihood (Positive Scale)

Calculates the positive log-likelihood contribution for movement process
error deviations under multiple IID structural assumptions. Deviations
are penalized as \\N(0, \sigma^2)\\ where \\\sigma\\ is drawn from
`PE_pars` according to the selected model structure. Under unstructured
movement, only origin-destination pairs that are adjacent (non-zero in
`adjacency_collapsed`) contribute to the likelihood; under CTMC movement
the deviations sit on each region's preference, and every region the map
keeps contributes.

## Usage

``` r
Get_move_PE_loglik(
  cont_vary_movement,
  PE_pars,
  move_devs,
  map_move_devs,
  do_recruits_move,
  adjacency_collapsed,
  move_type
)
```

## Arguments

- cont_vary_movement:

  Character string specifying the movement process error structure,
  `"iid_"` followed by the dims the deviations vary over, any of p
  (population), y (year), seas (season), a (age), s (sex): `"iid_y"` is
  one \\\sigma\\ per origin region, `"iid_p_y_seas_a_s"` one per
  population, origin region, season, age and sex. A dim left out shares
  one deviation, and one \\\sigma\\, across it.

- PE_pars:

  Array of movement process error parameters (log standard deviations)
  dimensioned `[pop, from_region, seas, age, sex]`. Exponentiated
  internally to obtain \\\sigma\\. Which dimensions are active depends
  on `cont_vary_movement`; unused dimensions should be fixed at a
  constant (e.g., index 1) via the parameter map.

- move_devs:

  Movement deviation array dimensioned
  `[pop, from_region, to_region, year, seas, age, sex]`. Under CTMC
  movement a deviation sits on a region's preference rather than on a
  pair, so `from_region` is that region and `to_region` has length one.

- map_move_devs:

  Integer array dimensioned
  `[pop, from_region, to_region, year, seas, age, sex]` mapping
  deviations to unique estimated parameters. Shared deviations hold the
  same integer value; dimensions are extracted from this array to
  determine loop bounds.

- do_recruits_move:

  Integer (0/1). If `0` and `n_ages >= 2`, age-1 recruits are excluded
  from the likelihood (loop starts at age 2). If `1`, all ages including
  recruits are penalized.

- adjacency_collapsed:

  `[n_regions x (n_regions - 1)]` matrix of allowable movement
  connections among regions, with self-retention collapsed out.
  Origin-destination pairs with a value of 0 are skipped and contribute
  nothing to the likelihood. Read only when `move_type == 0`.

- move_type:

  Integer specifying the movement formulation:

  - **0** = Unstructured multinomial logit movement

  - **1** = CTMC-based movement

  Decides whether `adjacency_collapsed` is read: an unstructured
  deviation belongs to a region pair, a CTMC deviation to a single
  region.

## Value

Numeric scalar: the positive log-likelihood contribution from movement
process error deviations. Negated externally to form the negative
log-likelihood.

## Details

**Note:** The returned value is on the *positive* log-likelihood scale.
It must be negated externally to form the negative log-likelihood.
