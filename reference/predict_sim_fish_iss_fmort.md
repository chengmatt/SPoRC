# Predict ISS fishery compositions under fishing mortality

Uses historical ISS fishery compositions and fishing mortality rates to
estimate ISS compositions in the projection year. Compositions are
scaled relative to the historical maximum fishing mortality with linear
interpolation between the minimum and maximum observed ISS values. If
historical values are not available, defaults to the mean or zero.

## Usage

``` r
predict_sim_fish_iss_fmort(ISS_FishComps, Fmort, y, sim)
```

## Arguments

- ISS_FishComps:

  Array of ISS fishery compositions with dimensions \`\[region, year,
  sex, fleet, sim\]\`.

- Fmort:

  Array of fishing mortality rates with dimensions \`\[region, year,
  fleet, sim\]\`.

- y:

  Integer, projection year index for prediction.

- sim:

  Integer, simulation index.

## Value

Array with predicted ISS values for year \`y\`.
