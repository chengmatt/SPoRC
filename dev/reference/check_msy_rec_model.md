# Reject MSY reference points for a mean recruitment fit

MSY is the maximum of equilibrium yield over a stock-recruit curve, so
it is only defined when the fit estimated one. A model fitted with
`rec_model = "mean_rec"` has `rec_model == 0`, and the equilibrium
recruitment helpers in `refpts_msy.R` branch on Ricker against
everything else, so an unguarded mean recruitment fit would be handed
Beverton-Holt reference points. Steepness is also mapped off under mean
recruitment, so those reference points would sit at the default
`h_trans = 0.6` rather than anything the model estimated. Both failures
are silent, which is why this is an error and not a warning.

## Usage

``` r
check_msy_rec_model(what, rec_model, sr_penalty = 0)
```

## Arguments

- what:

  Character. The requested reference point, after the deprecated `BH_`
  names have been mapped to their current equivalents.

- rec_model:

  Integer. `0` mean recruitment, `1` Beverton-Holt, `2` Ricker.

- sr_penalty:

  Integer. `0` none, `1` Beverton-Holt, `2` Ricker. Only meaningful
  under `rec_model == 0`, where it adds a curve penalized against the
  recruitment deviations.

## Value

`invisible(NULL)`. Called for the error.

## Details

SPR reference points do not go through the curve at all. They scale
spawning biomass per recruit by mean recruitment taken from `rep$Rec`,
so they stay well defined under every `rec_model` and are the right
request here.
