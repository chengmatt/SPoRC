# Map survey catchability parameters

Internal helper called by
[`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md)
to construct the TMB/RTMB factor map for `ln_srv_q`
`[n_regions × max_q_blocks × n_srv_fleets]`, the log-scale survey
catchability parameters. Parameters are only activated for region–fleet
combinations where survey index data are used
(`sum(UseSrvIdx[r,,,f]) > 0`)
`sum(input_list$data$UseSrvIdx_pop[,r,,,f] > 0)`; unused combinations
are mapped to `NA`.

## Usage

``` r
do_srv_q_mapping(input_list, srv_q_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$srv_q_blocks`,
  `$data$UseSrvIdx`, and `$data$UseSrvIdx_pop`.

- srv_q_spec:

  Character vector `[n_srv_fleets]` specifying the sharing structure for
  catchability. One of:

  `"est_all"`

  :   Separate catchability per region and block.

  `"est_shared_r"`

  :   Single catchability per block shared across all regions. Block
      membership is checked per region before assignment.

  `"fix"`

  :   All catchability parameters fixed at starting values (mapped to
      `NA`).

## Value

The input `input_list` with `$map$ln_srv_q` set to a factor vector.
Active parameters receive sequential integer indices; unused
region–fleet combinations and fixed parameters are `NA`.
