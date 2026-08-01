# Beta prior on tag reporting rate

Called once from the "Tag Reporting Rate (Prior)" section of
`SPoRC_rtmb.R`. Supports both a symmetric-beta parameterization
(`type == 0`) and a mean/sd beta parameterization (`type == 1`).

## Usage

``` r
get_tagrep_prior(conv_tag_fishrep_prior, conv_tag_fish_reporting_pars)
```

## Arguments

- conv_tag_fishrep_prior:

  Data frame with columns `region`, `block`, `fleet`, `type` (0 =
  symmetric beta, 1 = mean/sd beta), `mu`, `sd`, one row per penalized
  parameter.

- conv_tag_fish_reporting_pars:

  Array `[region, block, fish_fleet]` of tag reporting rate parameters
  on the logit scale.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `conv_tag_fishrep_prior`.
