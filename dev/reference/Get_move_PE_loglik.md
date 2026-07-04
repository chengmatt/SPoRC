# Compute Movement Process Error Log-Likelihood (Positive Scale)

Calculates the positive log-likelihood contribution for movement process
error deviations under multiple IID structural assumptions. Deviations
are penalized as \\N(0, \sigma^2)\\ where \\\sigma\\ is drawn from
`PE_pars` according to the selected model structure. Only
origin-destination pairs that are adjacent (non-zero in
`adjacency_collapsed`) contribute to the likelihood.

## Usage

``` r
Get_move_PE_loglik(
  PE_model,
  PE_pars,
  move_devs,
  map_move_devs,
  do_recruits_move,
  adjacency_collapsed,
  move_type
)
```

## Arguments

- PE_model:

  Integer specifying the movement process error structure. All models
  are IID; they differ in which dimensions share a common standard
  deviation. Models 1–5 are single-population (fix `pop = 1`); models
  6–10 estimate separate parameters per population:

  - **1** – IID across years (single \\\sigma\\ per origin region)

  - **2** – IID across ages (single \\\sigma\\ per origin region and
    age)

  - **3** – IID across years and ages

  - **4** – IID across years, ages, and sexes

  - **5** – IID across years, seasons, ages, and sexes

  - **6** – IID across populations and years

  - **7** – IID across populations and ages

  - **8** – IID across populations, years, and ages

  - **9** – IID across populations, years, ages, and sexes

  - **10** – IID across populations, years, seasons, ages, and sexes

- PE_pars:

  Array of movement process error parameters (log standard deviations)
  dimensioned `[pop, from_region, seas, age, sex]`. Exponentiated
  internally to obtain \\\sigma\\. Which dimensions are active depends
  on `PE_model`; unused dimensions should be fixed at a constant (e.g.,
  index 1) via the parameter map.

- move_devs:

  Movement deviation array dimensioned
  `[pop, from_region, to_region, year, seas, age, sex]`.

- map_move_devs:

  Integer array dimensioned
  `[pop, from_region, to_region, year, seas, age, sex]` mapping
  deviations to unique estimated parameters. Shared deviations carry the
  same integer value; dimensions are extracted from this array to
  determine loop bounds.

- do_recruits_move:

  Integer (0/1). If `0` and `n_ages >= 2`, age-1 recruits are excluded
  from the likelihood (loop starts at age 2). If `1`, all ages including
  recruits are penalized.

- adjacency_collapsed:

  Square `[n_regions x n_regions]` matrix of allowable movement
  connections among regions, excluding self-retention (diagonal entries
  should be 0). Origin-destination pairs with a value of 0 are skipped
  and contribute nothing to the likelihood.

- move_type:

  Integer specifying the movement formulation:

  - **0** = Unstructured multinomial logit movement

  - **1** = CTMC-based movement

  Currently used for dispatch context; likelihood computation is
  identical across movement types within this function.

## Value

Numeric scalar: the positive log-likelihood contribution from movement
process error deviations. Negated externally to form the negative
log-likelihood.

## Details

**Note:** The returned value is on the *positive* log-likelihood scale.
It must be negated externally to form the negative log-likelihood.
