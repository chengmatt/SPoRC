# Map fishery tag reporting rate parameters

Internal helper called by
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
to construct the TMB/RTMB factor map for `conv_tag_fish_reporting_pars`
`[n_regions × max_tagrep_blocks × n_fish_fleets]`, the logit-scale tag
reporting rate parameters. Parameters are only active for fleets with
`use_conv_fish_tagging = 1` and likelihoods that use reporting rates
(Poisson, negative binomial, or Dirichlet-multinomial conditioned on
releases: codes `0`, `1`, `2`, `4`).

## Usage

``` r
do_conv_tag_fish_reporting_pars_mapping(input_list, conv_tagrep_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$use_conv_fish_tagging`, `$data$conv_fish_tag_like`,
  `$data$n_regions`, `$data$n_fish_fleets`, and
  `$data$conv_tag_fish_reporting_blocks`.

- conv_tagrep_spec:

  Character string specifying the sharing structure for reporting rate
  parameters. One of:

  `"est_all"`

  :   Separate parameter per region, fleet, and block.

  `"est_shared_r"`

  :   Shared across regions within each fleet and block. Requires
      identical block structure across regions within each fleet.

  `"est_shared_f"`

  :   Shared across active fleets within each region and block. Requires
      identical block structure across fleets within each region.

  `"est_shared_r_f"`

  :   Shared across all regions and active fleets within each block.
      Requires identical block structure across all regions and fleets.

  `"fix"`

  :   All reporting rate parameters fixed at starting values (mapped to
      `NA`).

## Value

The input `input_list` with `$map$conv_tag_fish_reporting_pars` set to a
factor vector of length `prod(dim(par$conv_tag_fish_reporting_pars))`.
Active parameters receive sequential integer indices; inactive fleet
slots and fixed parameters are `NA`.

## Details

Sharing options that impose region- or fleet-wide pooling
(`"est_shared_r"`, `"est_shared_f"`, `"est_shared_r_f"`) require that
all pooled cells share an identical block structure; a mismatch raises
an error. When all fleets have tagging disabled, all parameters are
mapped to `NA`.
