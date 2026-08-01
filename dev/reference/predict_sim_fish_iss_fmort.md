# Predict fishery and discarded ISS under projected fishing mortality

Scales fishery input sample sizes for the projection year `y` based on
the relationship between fishing mortality and historical ISS values.
For each region-sex-fleet cell, the minimum and maximum ISS from the
conditioning period (`1:(y-1)`) are identified from years with positive,
non-NA values, and the projected ISS is obtained by linear interpolation
between those bounds using the ratio of projected \\F_y\\ to the
historical maximum \\F\\ (capped at 1). If no valid historical
observations exist for a cell, ISS is set to zero. If conditions for
scaling are not met (e.g., maximum historical \\F = 0\\), the mean
historical ISS is used as a fallback. All prior years (`1:(y-1)`) are
carried over unchanged from `ISS_FishComps`.

## Usage

``` r
predict_sim_fish_iss_fmort(ISS_FishComps, Fmort, y, seas, sim)
```

## Arguments

- ISS_FishComps:

  Array of fishery ISS values
  `[n_regions × n_yrs × n_seas × n_sexes × n_fish_fleets × n_sims]`.

- Fmort:

  Array of fishing mortality rates
  `[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]`.

- y:

  Integer. Projection year index for which ISS is predicted.

- seas:

  Integer. Season index.

- sim:

  Integer. Simulation replicate index.

## Value

Array `[n_regions × y × 1 × n_sexes × n_fish_fleets]` with historical
ISS values filled in for years `1:(y-1)` and the predicted ISS in year
`y`.
